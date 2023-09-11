import sys
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

    vcf = vcf[0]
    
    if len(vcf) == 0:
        raise ValueError("At least one VCF file must be supplied")
    if prefix is None and len(vcf) > 1:
        raise ValueError("If there is more than one vcf file, --prefix must be set")
    check_missing_file(vcf, die=True)

    prefix = create_prefix(vcf[0], prefix)

    # write annotations to excel file
    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
    logging.info("saving xlsx report: " + xlsx_name)
    writer = pandas.ExcelWriter(xlsx_name, engine='xlsxwriter')
    annotated_variants1.to_excel(writer, sheet_name='Variants', index=False)
    writer.close()
    
    csv_name = os.path.join(out_folder_path, prefix + ".annotated_variants.csv")
    logging.info("saving csv report: " + csv_name)
    annotated_variants1.to_csv(csv_name, index=False, )