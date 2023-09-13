import sys
import pysam
import logging
import gzip
import pandas
import os.path
import xlsxwriter

# LOGGING
logger = logging.getLogger(__name__)


def do_report(debug, vcf, prefix=None, min_vaf=0.0, out_folder_path="."):
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

    excel_headers = [
        "SAMPLE",
        "CHR",
        "POS",
        "REF",
        "ALT",
        "QUAL",
        "FILTER",
        "MQM_INFO",
        "MQMR_INFO",
        "QA_INFO",
        "QR_INFO",
        "SAF_INFO",
        "SAR_INFO",
        "SRF_INFO",
        "SRR_INFO",
        "SBR_INFO",
        "SBA_INFO",
        "POS_FILTER",
        "SBR_FILTER",
        "SBA_FILTER",
        "MQMR_FILTER",
        "AQR_FILTER",
        "GT_FORMAT",
        "QR_FORMAT",
        "AQR_FORMAT",
        "QA_FORMAT",
        "AQA_FORMAT",
        "INFO",
        "FORMAT",
    ]

    excel_table = [excel_headers]

    info_column_names = [
        "MQM",
        "MQMR",
        "QA",
        "QR",
        "SAF",
        "SAR",
        "SRF",
        "SRR",
        "SBR",
        "SBA",
    ]

    filter_column_names = ["POS", "SBR", "SBA", "MQMR", "AQR"]
    format_column_names = ["GT", "QR", "AQR", "QA", "AQA"]

    annotated_vcf = pysam.VariantFile(vcf)
    for variant in annotated_vcf.fetch():
        # complete INFO field string
        info_field_array = []
        for key, value in variant.info.items():
            if isinstance(value, tuple):
                value = value[0]
            info_field_array.append(f"{key}={value}")
        info_field_string = ";".join(info_field_array)

        # samples
        for sample in variant.samples.values():
            # skip sample if the VAF is too low
            if float(sample["VAF"][0]) <= float(min_vaf):
                continue

            # complete FORMAT field
            format_field_array = []
            for value in sample.values():
                if isinstance(value, tuple):
                    value = value[0]
                format_field_array.append(str(value))
            format_field_string = ":".join(format_field_array)

            row = []
            row.append(sample.name)

            # VCF
            vcf_start_columns = [
                variant.chrom,
                variant.pos,
                variant.ref,
                variant.alts,
                variant.qual,
                variant.filter,
            ]
            row.extend(vcf_start_columns)

            # INFO field
            info_column_data = []
            for field in info_column_names:
                info_column_data.append(variant.info[field])

            row.extend(info_column_data)

            # FORMAT filters
            filter_column_data = []
            for filter_name in filter_column_names:
                filter_column_data.append(
                    "PASS" if sample[filter_name] == 1 else "FAIL"
                )

            row.extend(filter_column_data)

            # FORMAT field
            format_column_data = []
            for format_name in format_column_names:
                format_column_data.append(sample[format_name])

            row.extend(format_column_data)

            # add complete INFO and FORMAT fields
            row.append(info_field_string)
            row.append(format_field_string)

            excel_table.append(row)

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    excel_pandasdf = pandas.read_table(excel_table)

    xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
    with pandas.ExcelWriter(xlsx_name, engine="xlsxwriter") as writer:
        excel_pandasdf.to_excel(writer, sheet_name="Variants", index=False)


if __name__ == "__main__":
    do_report(
        False,
        "tests/ashkenazim-trio/output2/vcfanno_annotated.vcf",
        None,
        0.01,
        "tests/random",
    )
