"""
Contains utility functions for mity modules.
"""
import logging
import os
import subprocess
import sys
from glob import glob
import pysam


class MityUtil:
    """
    Contains utility functions for mity modules.
    """

    MITY_DIR = "mitylib"
    GENOME_FILE = "mitylib/reference/b37d5.genome"
    REF_DIR = "reference"
    ANNOT_DIR = "annot"

    @classmethod
    def get_mity_dir(cls):
        """
        Get the directory path of the Mity library.

        Returns:
            str: The path to the Mity library directory.
        """
        path = os.path.dirname(sys.modules["mitylib"].__file__)
        return path

    @classmethod
    def tabix(cls, f: str):
        """
        Generate a tabix index for a bgzipped file.

        Parameters:
            f (str): The path to a bgzip compressed file.

        Returns:
            None
        """
        tabix_call = "tabix -f " + f
        logging.debug(tabix_call)
        subprocess.run(tabix_call, shell=True, check=False)

    @classmethod
    def select_reference_fasta(cls, reference: str, custom_reference_fa: str = None):
        """
        Select the reference genome fasta file.

        Parameters:
            reference (str): One of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10.
            custom_reference_fa (str, optional): The path to a custom reference genome, or None.

        Returns:
            str: The path to the selected reference genome fasta file.
        """
        if custom_reference_fa is not None and os.path.exists(custom_reference_fa):
            res = custom_reference_fa
        else:
            ref_dir = os.path.join(cls.get_mity_dir(), cls.REF_DIR)
            res = glob(f"{ref_dir}/{reference}.*.fa")
            logging.debug(",".join(res))
            assert len(res) == 1
            res = res[0]
        return res

    @classmethod
    def select_reference_genome(
        cls, reference: str, custom_reference_genome: str = None
    ):
        """
        Select the reference genome .genome file.

        Parameters:
            reference (str): One of the inbuilt reference genomes: hs37d5, hg19, hg38, mm10.
            custom_reference_genome (str, optional): The path to a custom reference .genome file, or None.

        Returns:
            str: The path to the selected reference .genome file.
        """
        if custom_reference_genome is not None and os.path.exists(
            custom_reference_genome
        ):
            res = custom_reference_genome
        else:
            ref_dir = os.path.join(cls.get_mity_dir(), cls.REF_DIR)
            logging.debug("Looking for .genome file in %s", ref_dir)
            res = glob(f"{ref_dir}/{reference}.genome")
            logging.debug(",".join(res))
            assert len(res) == 1
            res = res[0]
        return res

    @classmethod
    def vcf_get_mt_contig(cls, vcf: str):
        """
        Get the mitochondrial contig name and length from a VCF file.

        Parameters:
            vcf (str): Path to a VCF file.

        Returns:
            tuple: A tuple of contig name as a str and length as an int.
        """
        r = pysam.VariantFile(vcf, "r")
        chroms = r.header.contigs
        mito_contig = set(["MT", "chrM"]).intersection(chroms)
        assert len(mito_contig) == 1
        mito_contig = "".join(mito_contig)
        return r.header.contigs[mito_contig].name, r.header.contigs[mito_contig].length

    @classmethod
    def get_annot_file(cls, annotation_file_path: str):
        """
        Get the path to an annotation file.

        Parameters:
            annotation_file_path (str): The name of the annotation file.

        Returns:
            str: The path to the annotation file.
        """
        mitylibdir = cls.get_mity_dir()
        path = os.path.join(mitylibdir, cls.ANNOT_DIR, annotation_file_path)
        assert os.path.exists(path)
        return path

    @classmethod
    def make_prefix(cls, vcf_path: str):
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

    @classmethod
    def gsort(cls, input_path: str, output_path: str, genome: str):
        """
        Run gsort.
        """
        gsort_cmd = f"gsort {input_path} {genome} | bgzip -cf > {output_path}"
        subprocess.run(gsort_cmd, shell=True, check=False)
        MityUtil.tabix(output_path)
