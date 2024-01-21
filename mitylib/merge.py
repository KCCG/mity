"""
MITY Merge

Combines a nuclear VCF with MITY by adding MITY MT variants.
"""

import logging
import os.path
import sys
import gzip
import pysam
import pysam.bcftools
from mitylib.util import MityUtil


logger = logging.getLogger(__name__)


class Merge:
    """
    Combines a nuclear VCF with a MITY VCF
    """

    def __init__(
        self,
        debug,
        nuclear_vcf_path,
        mity_vcf_path,
        genome,
        output_dir=".",
        prefix=None,
        keep=False,
    ):
        self.debug = debug
        self.nuclear_vcf_path = nuclear_vcf_path
        self.mity_vcf_path = mity_vcf_path
        self.genome = genome

        self.bcftools_isec_path = ""
        self.bcftools_concat_path = ""
        self.merged_sorted_vcf_path = ""
        self.merged_unsorted_vcf_path = ""

        self.output_dir = output_dir
        self.prefix = prefix
        self.keep = keep

        self.run()

    def run(self):
        """
        Run mity merge.
        """

        if self.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Entered debug mode.")
        else:
            logger.setLevel(logging.INFO)

        self.run_checks()
        self.set_strings()
        self.run_bcftools_isec()
        self.run_bcftools_concat()
        self.write_merged()
        MityUtil.gsort(
            self.merged_unsorted_vcf_path, self.merged_sorted_vcf_path, self.genome
        )
        self.remove_intermediate_files()

    def run_checks(self):
        """
        Check whether the nuclear and mity vcfs are compatible based on:
            - contig matching
        """

        mity_contig = MityUtil.vcf_get_mt_contig(self.mity_vcf_path)
        nuclear_contig = MityUtil.vcf_get_mt_contig(self.nuclear_vcf_path)

        if mity_contig != nuclear_contig:
            logger.error("The VCF files use mitochondrial contigs.")
            sys.exit()

    def set_strings(self):
        """
        Sets:
            - bcftools_isec_path
            - bcftools_concat_path
            - merged_unsorted_vcf_path
            - merged_sorted_vcf_path
        """
        if self.prefix is None:
            self.prefix = MityUtil.make_prefix(self.mity_vcf_path)

        self.bcftools_isec_path = os.path.join(self.output_dir, "0000.vcf")
        self.bcftools_concat_path = os.path.join(
            self.output_dir, self.prefix + ".bcftools.concat.vcf"
        )

        self.merged_unsorted_vcf_path = os.path.join(
            self.output_dir, self.prefix + ".mity.merge.unsorted.vcf"
        )
        self.merged_sorted_vcf_path = os.path.join(
            self.output_dir, self.prefix + ".mity.merge.vcf.gz"
        )

    def run_bcftools_isec(self):
        """
        Run bcftools isec
        """

        with open(self.bcftools_isec_path, "w", encoding="utf-8") as file:
            print(
                pysam.bcftools.isec(
                    "-p",
                    self.output_dir,
                    "-C",
                    self.nuclear_vcf_path,
                    self.mity_vcf_path,
                ),
                end="",
                file=file,
            )

    def run_bcftools_concat(self):
        """
        Run bcftools concat
        """

        with open(self.bcftools_concat_path, "w", encoding="utf-8") as file:
            print(
                pysam.bcftools.concat(self.bcftools_isec_path, self.mity_vcf_path),
                end="",
                file=file,
            )

    def merge_description(self, nuclear_description: str, mity_description: str):
        """
        Merge nuclear and mity header description.
        """
        return f"If CHR=MT OR CHR=chrM: {mity_description}, otherwise: {nuclear_description}"

    def get_header_line_nums(self):
        """
        Returns a dictionary with key value pairs of ID : current position.

        E.g.
        ##FORMAT=<ID=AD,Number=.,Type=Integer,Description="Allelic depths ...">
        ##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Approximate rea...">
        ##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">
        ##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">

        ##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in g...">
        ##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency, f...">
        ##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of a...">
        ##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description="Z-score f...">


        Gives the dictionary:
        dict = {
            FORMAT: {
                AD: 0,      # this represents the line starting from 0 bytes
                DP: 80,     # this line starts from 80 bytes
                GQ: 160,
                GT: ...
            },
            INFO: {
                AC: ...,
                AF: ...,
                AN: ...,
                BaseQRankSum: ...
            }
        }

        This is then used to file.seek later. Note that this is a workaround to
        reduce memory footprint (it only stores position), but also give
        O(n + m) performance instead of O(nm) if we simply loop over each file
        from the beginning each time.
        """

        line_nums = {"FORMAT": {}, "INFO": {}}

        with gzip.open(self.mity_vcf_path, "rt") as file:
            position = 0
            file.seek(0)

            while True:
                position = file.tell()

                line = file.readline()

                if not line:
                    break

                if not line.startswith("##"):
                    break

                if line.startswith("##FORMAT") or line.startswith("##INFO"):
                    section, field_id = self.get_header_line_info(line)
                    line_nums[section][field_id] = position

        return line_nums

    def get_header_line_info(self, line: str):
        """
        Get section (INFO/FORMAT) and ID of a header line.
        """
        if line.startswith("##FORMAT"):
            section = "FORMAT"
        else:
            section = "INFO"

        line = line.split(",")[0]
        field_id = line.split("=")[-1]

        return section, field_id

    def get_header_description(self, line: str):
        """
        Get header description.
        """
        if "Description=" not in line:
            return ""

        return line.split("Description=")[1].strip('>" ')

    def make_new_line(self, nuclear_line, mity_file_line_num):
        """
        Makes an updated line with merged description.
        """
        with gzip.open(self.mity_vcf_path, "rt") as file:
            file.seek(mity_file_line_num)
            mity_line = file.readline()

        mity_line = mity_line.strip("\n")
        nuclear_line = nuclear_line.strip("\n")

        nuclear_description = self.get_header_description(nuclear_line)
        mity_description = self.get_header_description(mity_line)

        new_description = self.merge_description(nuclear_description, mity_description)

        new_line = nuclear_line.split("Description=")[0]
        new_line += f'Description="{new_description}">\n'

        return new_line

    def write_merged(self):
        """
        Make a new file with updated headers and add variants.
        """

        header_dict = self.get_header_line_nums()

        with open(self.merged_unsorted_vcf_path, "w", encoding="utf-8") as merged_file:
            with open(self.bcftools_concat_path, "r") as concat_file:
                for line in concat_file:
                    if line.startswith("##FORMAT") or line.startswith("##INFO"):
                        section, field_id = self.get_header_line_info(line)
                        if field_id in header_dict[section]:
                            new_line = self.make_new_line(
                                line, header_dict[section][field_id]
                            )
                            merged_file.write(new_line)
                        else:
                            merged_file.write(line)
                    else:
                        merged_file.write(line)

    def remove_intermediate_files(self):
        """
        Removes extra files from bcftools isec:
            README.txt
            sites.txt

        Optionally removes:
            bcftools.isec file
            bcftools.concat file
        """
        os.remove(os.path.join(self.output_dir, "README.txt"))
        os.remove(os.path.join(self.output_dir, "sites.txt"))

        if not self.keep:
            os.remove(self.merged_unsorted_vcf_path)
            os.remove(self.bcftools_isec_path)
            os.remove(self.bcftools_concat_path)
