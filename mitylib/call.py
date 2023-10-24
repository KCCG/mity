"""Mitochondrial variant calling."""
import subprocess
import logging
import os.path
import sys

import pysam

from mitylib.normalise import Normalise
from mitylib.util import MityUtil

logger = logging.getLogger(__name__)


class Call:
    """
    Mity call.
    """

    MIN_MQ = 30
    MIN_BQ = 24
    MIN_AF = 0.01
    MIN_AC = 4
    P_VAL = 0.002

    def __init__(
        self,
        debug,
        files,
        reference,
        genome=None,
        prefix=None,
        min_mq=MIN_MQ,
        min_bq=MIN_BQ,
        min_af=MIN_AF,
        min_ac=MIN_AC,
        p=P_VAL,
        normalise=True,
        out_folder_path=".",
        region=None,
    ):
        self.debug = debug
        self.files = files[0]
        self.reference = reference
        self.genome = genome
        self.prefix = prefix
        self.min_mq = min_mq
        self.min_bq = min_bq
        self.min_af = min_af
        self.min_ac = min_ac
        self.p = p
        self.normalise = normalise
        self.out_folder_path = out_folder_path
        self.region = region

        self.file_string = ""
        self.normalise_output_file = ""
        self.call_output_file = ""

        self.mity_cmd = ""
        self.sed_cmd = ""

        self.run()

    def run(self):
        """
        Run mity call.
        """

        if self.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Entered debug mode.")
        else:
            logger.setLevel(logging.INFO)

        self.run_checks()
        self.set_strings()
        self.set_region()
        self.set_mity_cmd()

        self.run_freebayes()

        if self.normalise:
            self.run_normalise()
        else:
            MityUtil.tabix(self.call_output_file)

    def run_normalise(self):
        """
        Run mity normalise.
        """
        logger.debug("Normalising and Filtering variants")

        try:
            Normalise(
                debug=self.debug,
                vcf=self.call_output_file,
                reference_fasta=self.reference,
                output_file=self.normalise_output_file,
                allsamples=False,
                p=self.p,
                genome=self.genome,
            )
        finally:
            os.remove(self.call_output_file)

    def run_freebayes(self):
        """
        Run freebayes.
        """
        freebayes_call = (
            f"set -o pipefail && freebayes -f {self.reference} {self.file_string} "
            f"--min-mapping-quality {self.min_mq} "
            f"--min-base-quality {self.min_bq} "
            f"--min-alternate-fraction {self.min_af} "
            f"--min-alternate-count {self.min_ac} "
            f"--ploidy 2 "
            f"--region {self.region} "
            f"| sed 's/##source/##freebayesSource/' "
            f"| sed 's/##commandline/##freebayesCommandline/' "
            f"| {self.sed_cmd} | bgzip > {self.call_output_file}"
        )

        logger.info("Running FreeBayes in sensitive mode")
        logger.debug(freebayes_call)
        res = subprocess.run(
            freebayes_call,
            shell=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            executable="/bin/bash",
            check=False,
        )
        logger.debug("Freebayes result code: %s", res.returncode)

        if res.returncode != 0:
            logger.error("FreeBayes failed: %s", res.stderr)
            exit(1)

        if os.path.isfile(self.call_output_file):
            logger.debug("Finished running FreeBayes")

    def set_mity_cmd(self):
        """
        Creates the mity command for embedding into the vcf.
        """

        mity_cmd = (
            f'##mityCommandline="mity call --reference {self.reference} --prefix {self.prefix} '
            f"--min-mapping-quality {self.min_mq} --min-base-quality {self.min_bq} "
            f"--min-alternate-fraction {self.min_af} --min-alternate-count {self.min_ac} "
            f"--out-folder-path {self.out_folder_path} --region {self.region}"
        )

        if self.normalise:
            mity_cmd += f" --normalise --p {self.p}"

        mity_cmd += " " + " ".join(self.files)
        mity_cmd += '"'
        mity_cmd = mity_cmd.replace("/", "\/")

        logger.debug(mity_cmd)

        # overwrite a redundant freebayes header line with the mity command line
        sed_cmd = f"sed 's/^##phasing=none/{mity_cmd}/g'"
        logger.debug(sed_cmd)

        self.mity_cmd = mity_cmd
        self.sed_cmd = sed_cmd

    def set_region(self):
        """
        Sets the region if not specified.
        """
        if self.region is None:
            self.region = self.bam_get_mt_contig(self.files[0], as_string=True)

    def set_strings(self):
        """
        Sets:
            prefix
            file_string
            call_output_file
            normalise_output_file
        """
        self.prefix = self.create_prefix(self.files[0], self.prefix)
        self.file_string = " ".join(["-b " + _file for _file in reversed(self.files)])
        self.normalise_output_file = os.path.join(
            self.out_folder_path, self.prefix + ".mity.normalise.vcf.gz"
        )
        self.call_output_file = os.path.join(
            self.out_folder_path, self.prefix + ".mity.call.vcf.gz"
        )

    def run_checks(self):
        """
        Check for valid input.
        """

        if len(self.files) > 1 and self.prefix is None:
            raise ValueError(
                "If there is more than one bam/cram file, --prefix must be set"
            )

        self.check_missing_file(self.files, die=True)
        self.prefix = self.create_prefix(self.files[0], self.prefix)

        if not all(map(self.bam_has_rg, self.files)):
            logger.error("At least one BAM/CRAM file lacks an @RG header")
            exit(1)

        if self.normalise and self.genome is None:
            logger.error("A genome file should be supplied if mity call normalize=True")
            sys.exit(1)

    def bam_has_rg(self, bam):
        """
        Does the BAM or CRAM File have an @RG header? This is critical for mity
        to correctly call variants.

        :param bam: str: path to bam or cram file
        :return: True/False
        >>> bam_has_RG('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam')
        """
        r = pysam.AlignmentFile(bam, "rb")
        return len(r.header["RG"]) > 0

    def bam_get_mt_contig(self, bam, as_string=False):
        """
        get the mitochondrial contig name and length from a BAM file
        :param bam: path to a bam or cram file
        :return: a tuple of contig name as str and length as int

        >>> bam_get_mt_contig('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam', False)
        ('chrM', 16569)
        >>> bam_get_mt_contig('NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam', True)
        'chrM:1-16569'
        """
        r = pysam.AlignmentFile(bam, "rb")
        chroms = [str(record.get("SN")) for record in r.header["SQ"]]
        mito_contig = {"MT", "chrM"}.intersection(chroms)
        assert len(mito_contig) == 1
        mito_contig = "".join(mito_contig)
        res = None
        for record in r.header["SQ"]:
            if mito_contig == record["SN"]:
                res = record["SN"], record["LN"]
        if res is not None and as_string:
            res = res[0] + ":1-" + str(res[1])
        return res

    def create_prefix(self, file_name, prefix=None):
        """
        Generate a prefix for Mity functions. If a custom prefix is not provided,
        the function uses the filename without the file extension (.vcf, .bam, .cram, .bed).

        Parameters:
            file_name (str): The filename, including extensions (e.g., .vcf, .bam, .cram, .bed).
            prefix (str, optional): An optional custom prefix. If None, the function generates a prefix from the file name.

        Returns:
            str: The generated or custom prefix for the Mity function.
        """
        supported_extensions = {".vcf": ".vcf", ".bam": ".bam", ".cram": ".cram"}

        ext = os.path.splitext(file_name)[1]
        if ext in supported_extensions:
            return (
                prefix
                if prefix is not None
                else os.path.basename(file_name).replace(ext, "")
            )
        else:
            raise ValueError("Unsupported file type")

    def check_missing_file(self, file_list, die=True):
        """
        Check if input files exist.
        """
        missing_files = []
        for item in file_list:
            if not os.path.isfile(item):
                missing_files.append(item)
        if die and len(missing_files) > 0:
            raise ValueError("Missing these files: " + ",".join(missing_files))
        return missing_files
