import sys
import pysam
import logging
import gzip
import pandas
import os.path
import xlsxwriter
from .util import check_missing_file, create_prefix, make_hgvs, get_annot_file

# LOGGING
logger = logging.getLogger(__name__)

def do_report(debug, vcf, prefix=None, min_vaf=0.0, out_folder_path = "."):
    """
    Create a mity report
    :param vcf: the path to a vcf file
    :param prefix: the optional prefix. This must be set if there is >1 vcf files
    :param min_vaf: only include vairants with vaf > min_vaf in the report
    :return:
    """
    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Entered debug mode.")
    else:
        logger.setLevel(logging.INFO)

    excel_table = []
    excel_headers = []
    annotated_vcf = pysam.VariantFile(vcf)
    for variant in annotated_vcf.fetch():
        for sample in variant.samples.values():
            
            if float(sample["VAF"]) > float(min_vaf):
                continue

            excel_row = []
            excel_row.append(sample.name)

            # FILTER
            

            # INFO field
            # FORMAT field


