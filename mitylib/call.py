"""Mitochondrial variant calling."""

import subprocess
import logging
import os.path
import sys
from typing import Any, Optional, Tuple
import urllib.request

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
        output_dir=".",
        region=None,
        bam_list=False,
        keep=False,
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
        self.output_dir = output_dir
        self.region = region
        self.bam_list = bam_list
        self.keep = keep

        self.file_string = ""
        self.normalised_vcf_path = ""
        self.call_vcf_path = ""

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

        if self.bam_list:
            self.get_files_from_list()
        self.run_checks()
        self.set_strings()
        self.set_region()
        self.set_mity_cmd()

        self.run_freebayes()

        if self.normalise:
            self.run_normalise()
        else:
            MityUtil.tabix(self.call_vcf_path)

    def run_normalise(self):
        """
        Run mity normalise.
        """
        logger.debug("Normalising and Filtering variants")

        try:
            Normalise(
                debug=self.debug,
                vcf=self.call_vcf_path,
                reference_fasta=self.reference,
                prefix=self.prefix,
                output_dir=self.output_dir,
                allsamples=False,
                p=self.p,
                genome=self.genome,
                keep=self.keep,
            )
        finally:
            if not self.keep:
                os.remove(self.call_vcf_path)

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
            f"| {self.sed_cmd} | bgzip > {self.call_vcf_path}"
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

        if os.path.isfile(self.call_vcf_path):
            logger.debug("Finished running FreeBayes")

    def set_mity_cmd(self):
        """
        Creates the mity command for embedding into the vcf.
        """

        mity_cmd = (
            f'##mityCommandline="mity call --reference {self.reference} --prefix {self.prefix} '
            f"--min-mapping-quality {self.min_mq} --min-base-quality {self.min_bq} "
            f"--min-alternate-fraction {self.min_af} --min-alternate-count {self.min_ac} "
            f"--out-folder-path {self.output_dir} --region {self.region}"
        )

        if self.normalise:
            mity_cmd += f" --normalise --p {self.p}"

        mity_cmd += " " + " ".join(self.files)
        mity_cmd += '"'
        mity_cmd = mity_cmd.replace("/", "\\/")

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
        if self.prefix is None:
            self.prefix = self.make_prefix(self.files[0], self.prefix)

        self.file_string = " ".join(["-b " + _file for _file in reversed(self.files)])

        self.normalised_vcf_path = os.path.join(
            self.output_dir, self.prefix + ".mity.normalise.vcf.gz"
        )
        self.call_vcf_path = os.path.join(
            self.output_dir, self.prefix + ".mity.call.vcf.gz"
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
        self.prefix = self.make_prefix(self.files[0], self.prefix)

        if not all(map(self.bam_has_rg, self.files)):
            logger.error("At least one BAM/CRAM file lacks an @RG header")
            exit(1)

        if self.normalise and self.genome is None:
            logger.error("A genome file should be supplied if mity call normalize=True")
            sys.exit(1)

    def bam_has_rg(self, bam):
        """
        Check whether a BAM or CRAM file contains a valid @RG header,
        which is critical for accurate variant calling with mity.

        Parameters:
            - bam (str): Path to a BAM or CRAM file.

        Returns:
            - bool: True if the file has a valid @RG header, False otherwise.
        """
        r = pysam.AlignmentFile(bam, "rb")
        return len(r.header["RG"]) > 0

    def bam_get_mt_contig(
        self, bam: str, as_string: bool = False
    ) -> Optional[Tuple[Any, Any] | str]:
        """
        Retrieve mitochondrial contig information from a BAM or CRAM file.

        Parameters:
            - bam (str): Path to a BAM or CRAM file.

        Keyword Arguments:
            - with_coordinates (bool): If True, the result includes the contig
            name with coordinates (e.g., 'chrM:1-16569'). If False, it returns a
            tuple containing the contig name (str) and its length (int).

        Returns:
            - If with_coordinates is False, a tuple (contig_name, contig_length).
            - If with_coordinates is True, a string with coordinates.
        """
        r = pysam.AlignmentFile(bam, "rb")
        chroms = [str(record.get("SN")) for record in r.header["SQ"]]
        mito_contig_intersection = {"MT", "chrM"}.intersection(chroms)

        assert len(mito_contig_intersection) == 1
        mito_contig = "".join(mito_contig_intersection)

        res = None

        for record in r.header["SQ"]:
            if mito_contig == record["SN"]:
                res = record["SN"], record["LN"]

        if res is not None and as_string:
            res = res[0] + ":1-" + str(res[1])

        return res

    def make_prefix(self, file_name: str, prefix: Optional[str] = None):
        """
        Generate a prefix for Mity functions if a custom prefix is not provided,
        the function uses the filename without the file extension (.vcf, .bam, .cram, .bed).

        Parameters:
            file_name (str): The filename, including extensions (e.g., .vcf, .bam, .cram, .bed).
            prefix (str, optional): An optional custom prefix. If None, the function generates a
                                    prefix from the file name.

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

    def get_files_from_list(self, die: bool = True):
        """
        Get the list of BAM / CRAM files from the provided file
        """
        if len(self.files) > 1:
            raise ValueError(
                "--bam-file-list Argument expects only 1 file to be provided."
            )
        with open(self.files[0], "r") as f:
            self.files = f.read().splitlines()

    def check_missing_file(self, file_list: str, die: bool = True):
        """
        Check if input files exist.
        """
        missing_files = []
        for item in file_list:
            if item.lower().startswith("http"):
                try:
                    urllib.request.urlopen(item).getcode()
                except:
                    missing_files.append(item)
            elif not os.path.isfile(item):
                missing_files.append(item)
        if die and len(missing_files) > 0:
            raise ValueError("Missing these files: " + ",".join(missing_files))
        return missing_files
