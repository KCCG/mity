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
    def tabix(cls, f):
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
    def select_reference_fasta(cls, reference, custom_reference_fa=None):
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
    def select_reference_genome(cls, reference, custom_reference_genome=None):
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
    def vcf_get_mt_contig(cls, vcf):
        """
        Get the mitochondrial contig name and length from a VCF file.

        Parameters:
            vcf (str): Path to a VCF file.

        Returns:
            tuple: A tuple of contig name as a str and length as an int.
        """
        r = pysam.VariantFile(vcf, "r", require_index=True)
        chroms = r.header.contigs
        mito_contig = set(["MT", "chrM"]).intersection(chroms)
        assert len(mito_contig) == 1
        mito_contig = "".join(mito_contig)
        return r.header.contigs[mito_contig].name, r.header.contigs[mito_contig].length

    @classmethod
    def get_annot_file(cls, f):
        """
        Get the path to an annotation file.

        Parameters:
            f (str): The name of the annotation file.

        Returns:
            str: The path to the annotation file.
        """
        mitylibdir = cls.get_mity_dir()
        p = os.path.join(mitylibdir, cls.ANNOT_DIR, f)
        assert os.path.exists(p)
        return p
