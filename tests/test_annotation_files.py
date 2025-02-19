import os
import glob
import pysam
import pytest

VCF_DIR = "/home/tzeng/tfiles/cci/projects/mity-project/mity/mitylib/annot_chrm/"


def get_vcf_files(directory):
    return glob.glob(os.path.join(directory, "*.vcf.gz"))


@pytest.mark.parametrize("vcf_file", get_vcf_files(VCF_DIR))
def test_vcf_file_validity(vcf_file):
    try:
        vcf = pysam.VariantFile(vcf_file)
        for record in vcf.fetch():
            assert record.info is not None
    except Exception as e:
        pytest.fail(f"Failed to parse {vcf_file}: {e}")
