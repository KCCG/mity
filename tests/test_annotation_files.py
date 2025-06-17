import os
import glob
import pysam
import pytest
import mitylib

ANNOT_CHRM_DIR = os.path.join(mitylib.__path__[0], "annot_chrm")
ANNOT_MT_DIR = os.path.join(mitylib.__path__[0], "annot_mt")


def get_vcf_files(directory):
    """
    Returns a list of all VCF files in a given directory.
    """
    return glob.glob(os.path.join(directory, "*.vcf.gz"))


@pytest.mark.parametrize("vcf_file", get_vcf_files(ANNOT_CHRM_DIR) + get_vcf_files(ANNOT_MT_DIR))
def test_vcf_file_validity(vcf_file: str) -> None:
    """
    Pytest test that verifies all VCF files in two directories are correctly formatted.
    """
    try:
        vcf = pysam.VariantFile(vcf_file)
        for record in vcf.fetch():
            assert record.info is not None
    except Exception as e:
        pytest.fail(f"Failed to parse {vcf_file}: {e}")


def get_bed_files(directory):
    """
    Returns a list of all BED files in a given directory.
    """
    return glob.glob(os.path.join(directory, "*.bed.gz"))


def validate_bed_file(bed_file: str) -> bool:
    """
    Validates a BED file by checking column count, numeric values, and formatting.

    Args:
        bed_file (str): Path to the BED file

    Returns:
        bool: True if valid, False otherwise
    """
    try:
        with pysam.TabixFile(bed_file, parser=pysam.asTuple()) as tbx:
            for row in tbx.fetch():
                cols = row[:3]  # At minimum, BED files must have 3 columns
                if len(cols) < 3:
                    return False

                chrom, start, end = cols

                if not chrom:
                    return False

                if not start.isdigit() or not end.isdigit():
                    return False

                start, end = int(start), int(end)

                if start < 0 or end < 0:
                    return False

                if end <= start:
                    return False

        return True
    except Exception:
        return False


@pytest.mark.parametrize("bed_file", get_bed_files(ANNOT_CHRM_DIR) + get_bed_files(ANNOT_MT_DIR))
def test_bed_file_format(bed_file):
    """
    Pytest test that verifies all BED files in two directories are correctly formatted.

    Args:
        bed_file (str): Path to the BED file (parameterized by pytest)
    """
    assert validate_bed_file(bed_file), f"Invalid BED file: {bed_file}"

def test_all_tbi_files_exist():
    """
    Pytest test that verifies all files in the annotation directories exist.
    """
    for directory in [ANNOT_CHRM_DIR, ANNOT_MT_DIR]:
        for file in get_vcf_files(directory) + get_bed_files(directory):
            assert os.path.exists(f"{file}.tbi"), f"tbi file does not exist for: {file}"