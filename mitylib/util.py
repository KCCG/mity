"""
Contains utility functions for mity modules.
"""

import logging
import os
import subprocess
import sys
from glob import glob
from typing import Optional, Tuple
import pysam


class MityUtil:
    """
    Contains utility functions for mity modules.
    """

    MITY_DIR = "mitylib"
    REF_DIR = "reference"
    ANNOT_DIR = "annot"

    @staticmethod
    def get_mity_dir():
        """
        Get the directory path of the Mity library.

        Returns:
            str: The path to the Mity library directory.
        """
        return os.path.dirname(sys.modules["mitylib"].__file__)

    @staticmethod
    def tabix(bgzipped_file: str) -> None:
        """
        Generate a tabix index for a bgzipped file.

        Parameters:
            bgzipped_file (str): The path to a bgzip compressed file.

        Returns:
            None
        """
        tabix_call = "tabix -f " + bgzipped_file
        logging.debug(tabix_call)
        subprocess.run(tabix_call, shell=True, check=False)

    @staticmethod
    def select_reference_fasta(
        reference: str, custom_reference_fa: Optional[str] = None
    ) -> str:
        """
        Select the reference genome fasta file.

        Parameters:
            reference (str): One of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10.
            custom_reference_fa (str, optional): The path to a custom reference genome, or None.

        Returns:
            str: The path to the selected reference genome fasta file.
        """
        if custom_reference_fa is not None:
            if not os.path.exists(custom_reference_fa):
                raise FileNotFoundError(
                    f"--custom-reference-fasta file: {custom_reference_fa} cannot be found."
                )
            return custom_reference_fa

        ref_dir = os.path.join(MityUtil.get_mity_dir(), MityUtil.REF_DIR)
        res = glob(f"{ref_dir}/{reference}.*.fa")
        logging.debug(",".join(res))
        assert len(res) == 1

        return res[0]

    @staticmethod
    def select_reference_genome(
        reference: str, custom_reference_genome: Optional[str] = None
    ) -> str:
        """
        Select the reference genome .genome file.

        Parameters:
            reference (str): One of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10.
            custom_reference_genome: The path to a custom reference .genome file, or None.

        Returns:
            str: The path to the selected reference .genome file.
        """
        if custom_reference_genome is not None:
            if not os.path.exists(custom_reference_genome):
                raise FileNotFoundError(
                    f"--custom-reference-genome file: {custom_reference_genome} cannot be found."
                )
            return custom_reference_genome

        ref_dir = os.path.join(MityUtil.get_mity_dir(), MityUtil.REF_DIR)
        logging.debug("Looking for .genome file in %s", ref_dir)
        res = glob(f"{ref_dir}/{reference}.genome")
        logging.debug(",".join(res))
        assert len(res) == 1
        return res[0]

    @staticmethod
    def vcf_get_mt_contig(vcf: str) -> Tuple[str, Optional[int]]:
        """
        Get the mitochondrial contig name and length from a VCF file.

        Parameters:
            vcf (str): Path to a VCF file.

        Returns:
            tuple: A tuple of contig name as a str and length as an int.
        """
        r = pysam.VariantFile(vcf, "r")
        chroms = r.header.contigs
        mito_contig_intersection = set(["MT", "chrM"]).intersection(chroms)

        assert len(mito_contig_intersection) == 1

        mito_contig = "".join(mito_contig_intersection)

        mt_contig_name = r.header.contigs[mito_contig].name
        mt_contig_length = r.header.contigs[mito_contig].length

        return (mt_contig_name, mt_contig_length)

    @staticmethod
    def get_annot_file(annotation_file_path: str):
        """
        Get the path to an annotation file.

        Parameters:
            annotation_file_path (str): The name of the annotation file.

        Returns:
            str: The path to the annotation file.
        """
        mitylibdir = MityUtil.get_mity_dir()
        path = os.path.join(mitylibdir, MityUtil.ANNOT_DIR, annotation_file_path)
        assert os.path.exists(path)
        return path

    @staticmethod
    def make_prefix(vcf_path: str):
        """
        Make a prefix based on the input vcf path. This handles vcf files from
        previous steps of mity. e.g. from call to normalise, etc.

        Format of MITY output filenames:
            prefix.mity.call.vcf.gz
            prefix.mity.normalise.vcf.gz
            prefix.mity.merge.vcf.gz
            prefix.report.xlsx
        """

        prefix = (
            os.path.basename(vcf_path)
            .replace(".mity", "")
            .replace(".call", "")
            .replace(".normalise", "")
            .replace(".merge", "")
            .replace(".report", "")
            .replace(".vcf.gz", "")
        )

        return prefix

    @staticmethod
    def gsort(input_path: str, output_path: str, genome: str):
        """
        Run gsort.
        """
        gsort_cmd = f"gsort {input_path} {genome} | bgzip -cf > {output_path}"
        subprocess.run(gsort_cmd, shell=True, check=False)
        MityUtil.tabix(output_path)
