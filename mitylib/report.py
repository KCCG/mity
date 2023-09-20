import sys
import pysam
import logging
import gzip
import pandas
import os.path
import xlsxwriter
import subprocess
from collections import namedtuple

# LOGGING
logger = logging.getLogger(__name__)

# def make_excel_row(variant, sample)


def make_hgvs(pos, ref, alt):
    """
    Creates HGVS syntax used in the HGVS column of the table/excel spreadsheet.
    """
    if len(alt) > 1 or len(ref) > 1:
        # this is an indel
        if len(ref) > len(alt):
            # this is a del
            delet = ref[1:]
            if len(delet) == 1:
                hgvs_pos = int(pos) + 1
            elif len(delet) > 1:
                hvgs_pos_start = int(pos) + 1
                hvgs_pos_end = int(pos) + len(delet)
                hgvs_pos = str(hvgs_pos_start) + "_" + str(hvgs_pos_end)
            hgvs = "m." + str(hgvs_pos) + "del"

        else:
            # this is an ins
            ins = alt[1:]
            if len(ins) == 1:
                hgvs_pos = int(pos) + 1
            elif len(ins) > 1:
                hvgs_pos_start = int(pos) + 1
                hvgs_pos_end = int(pos) + len(ins)
                hgvs_pos = str(hvgs_pos_start) + "_" + str(hvgs_pos_end)
            hgvs = "m." + str(hgvs_pos) + "ins"

    else:
        # this is a SNP
        hgvs = "m." + str(pos) + str(ref) + ">" + str(alt)
    return hgvs


def vcfanno_call(vcf):
    """
    Calls vcfanno to annotate the output of mity normalise with the relevant
    annotations based on vcfanno-config.toml in mitylib/
    """
    logger.debug("Running vcfanno...")

    # annotated_file name
    annotated_file = vcf.replace(".vcf.gz", ".mity.annotated.vcf")

    # vcfanno call
    vcfanno_cmd = f"vcfanno -p 4 vcfanno-config.toml {vcf} > {annotated_file}"
    res = subprocess.run(
        vcfanno_cmd,
        shell=True,
        check=False,
        # capture the output of the vcfanno command since it tends to produce
        # long warning messages, output shown in --debug mode
        capture_output=True,
        text=True,
    )

    logger.debug("vcfanno output:")
    logger.debug(res.stdout)

    return annotated_file


def make_info_string(variant):
    """
    Takes the INFO fields from the variant and recreates the format as it
    appears in a vcf file. i.e.

    FIELD=X;FIELD=Y;etc
    """
    info_field_array = []
    for key, value in variant.info.items():
        if isinstance(value, tuple):
            value = value[0]
        info_field_array.append(f"{key}={value}")
    info_field_string = ";".join(info_field_array)
    return info_field_string


def make_format_string(sample):
    """
    Takes the FORMAT fields from the sample and recreates the format as it
    appears in a vcf file. i.e.

    FIELD=X;FIELD=Y;etc
    """
    format_field_array = []
    for value in sample.values():
        if isinstance(value, tuple):
            value = value[0]
        format_field_array.append(str(value))
    format_field_string = ":".join(format_field_array)
    return format_field_string


def make_headers():
    """
    Returns a list that contains the headers of the excel spreadsheet in order
    """

    # separating the excel headers into sections to make it easier to read
    headers_dict = {
        "start": [
            "SAMPLE",
            "HGVS",
            "GENE/LOCUS",
            "GENE/LOCUS DESCRIPTION",
            "TOTAL LOCUS DEPTH",
            "VARIANT HETEROPLASMY",
            "ALT DEPTH",
            "REF DEPTH",
            "TOTAL SAMPLE DEPTH",
            "VARIANT QUALITY",
            "TIER",
            "COHORT COUNT",
            "COHORT FREQUENCY",
        ],
        "basic": [
            "CHR",
            "POS",
            "REF",
            "ALT",
            "QUAL",
            "FILTER",
        ],
        "info": [
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
        ],
        "filters": [
            "POS_FILTER",
            "SBR_FILTER",
            "SBA_FILTER",
            "MQMR_FILTER",
            "AQR_FILTER",
        ],
        "annotations": [
            "COMMERCIAL PANELS",
            "PHYLOTREE HAPLOTYPE",
            "MITOTIP SCORE",
            "MITOTIP PERCENTILE",
            "MITOTIP INTERPRETATION",
            "ANTICODON",
            "ALLELE FREQUENCY MITOMAP",
            "DISEASE MITOMAP",
            "MGRB FREQUENCY",
            "MGRB FILTER",
            "MGRB AC",
            "MGRB AN",
            "PHYLOTREE MUT",
            "LOCUS MITOMAP",
            "NUMBER OF REFERENCES MITOMAP",
            "VARIANT AMINO ACID CHANGE MITOMAP",
            "CODON POSITION MITOMAP",
            "CODON NUMBER MITOMAP",
            "NUM DISEASE REFERENCES MITOMAP",
            "RNA MITOMAP",
            "HOMOPLASMY MITOMAP",
            "HETEROPLASMY MITOMAP",
            "STATUS MITOMAP",
            "DISEASE AMINO ACID CHANGE MITOMAP",
            "GENBANK FREQUENCY MITOMAP",
        ],
        "format": [
            "GT_FORMAT",
            "QR_FORMAT",
            "AQR_FORMAT",
            "QA_FORMAT",
            "AQA_FORMAT",
        ],
        "vcf_info_format": [
            "INFO",
            "FORMAT",
        ],
    }

    excel_headers = []
    for header_group in headers_dict.values():
        excel_headers += header_group

    return excel_headers


def clean_string(s):
    if isinstance(s, tuple):
        s = s[0]

    s = str(s)
    s = s.replace('"', "").replace("(", "").replace(")", "")

    return s


def make_table(annotated_file, min_vaf):
    """
    Takes a vcfanno annotated vcf file and returns a formatted table with
    relevant information. The vcf names are hardcoded in vcf_column_names
    """

    vcf_column_names = {
        # Note that "start_info/format" in vcf_column_names does not include
        # some elements in headers_dict such as "SAMPLE" and "HGVS" since they
        # do not directly map to a field in the VCF, and are calculated
        # separately.
        "start_info": ["Map_Locus", "Description", "DP"],
        "start_format": ["VAF", "AO", "RO", "DP", "q", "tier"],
        "info": [
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
        ],
        "filters": ["POS", "SBR", "SBA", "MQMR", "AQR"],
        "format": ["GT", "QR", "AQR", "QA", "AQA"],
        "annotations": [
            "commercial_panels",
            "phylotree_haplotype",
            "MitoTip_score",
            "MitoTip_percentile",
            "MitoTip_interpretation",
            "anticodon",
        ],
        # Note that annotations are split into two parts because of "ALLELE
        # FREQUENCY MITOMAP which needs to be calculated separately"
        "annotations2": [
            "disease_mitomap",
            "MGRB_frequency",
            "MGRB_FILTER",
            "MGRB_AC",
            "MGRB_AN",
            "phylotree_mut",
            "locus_mitomap",
            "num_references_mitomap",
            "variant_amino_acid_change_mitomap",
            "codon_position_mitomap",
            "codon_number_mitomap",
            "num_disease_references_mitomap",
            "RNA_mitomap",
            "homoplasmy_mitomap",
            "heteroplasmy_mitomap",
            "status_mitomap",
            "disease_amino_acid_change_mitomap",
            "GenBank_frequency_mitomap",
        ],
    }

    excel_table = []

    # cohort count record keeping
    VariantDef = namedtuple("VariantDef", ["pos", "alt"])
    cohort_count = {}

    annotated_vcf = pysam.VariantFile(annotated_file)
    for variant in annotated_vcf.fetch():
        vcf_start_columns = [
            variant.chrom,
            variant.pos,
            variant.ref,
            variant.alts[0],
            variant.qual,
            variant.filter.keys()[0],
        ]

        # samples
        for sample in variant.samples.values():
            # skip sample if the VAF is too low
            if float(sample["VAF"][0]) <= float(min_vaf):
                continue

            # if the sample passes the VAF filter, add to the cohort count
            if VariantDef(variant.pos, variant.alts[0]) in cohort_count:
                cohort_count[VariantDef(variant.pos, variant.alts[0])] += 1
            else:
                cohort_count[VariantDef(variant.pos, variant.alts[0])] = 1

            row = []

            # headers_dict -> start
            row.append(sample.name)
            row.append(make_hgvs(variant.pos, variant.ref, variant.alts[0]))

            for name in vcf_column_names["start_info"]:
                if name in variant.info.keys():
                    row.append(clean_string(variant.info[name]))
                else:
                    row.append(".")

            for name in vcf_column_names["start_format"]:
                if name in sample.keys():
                    row.append(clean_string(sample[name]))
                else:
                    row.append(".")

            # placeholder for cohort count and cohort frequency
            row.append(0)
            row.append(0)

            # headers_dict -> basic
            row.extend(vcf_start_columns)

            # headers_dict -> info
            info_column_data = []
            for field in vcf_column_names["info"]:
                info_column_data.append(clean_string(variant.info[field]))

            row.extend(info_column_data)

            # headers_dict -> filters
            filter_column_data = []
            for filter_name in vcf_column_names["filters"]:
                filter_column_data.append(
                    # Note that the field is a tuple (1,) so we take the first
                    # index since there is only ever one value (i.e. true/false)
                    "PASS"
                    if sample[filter_name + "_filter"][0] == 1
                    else "FAIL"
                )

            row.extend(filter_column_data)

            # headers_dict -> annotations

            for name in vcf_column_names["annotations"]:
                if name in variant.info.keys():
                    row.append(clean_string(variant.info[name]))
                else:
                    row.append(".")

            # calculate ALLELE FREQUENCY MITOMAP
            # TODO: check this formula
            if "Genbank_frequency_mitomap" in sample.keys():
                allele_frequency_mitomap = round(
                    float(sample["GenBank_frequency_mitomap"]) / 32050, 3
                )
            else:
                allele_frequency_mitomap = "."

            row.append(allele_frequency_mitomap)

            for name in vcf_column_names["annotations2"]:
                if name in variant.info.keys():
                    row.append(clean_string(variant.info[name]))
                else:
                    row.append(".")

            # headers_dict -> format
            format_column_data = []
            for format_name in vcf_column_names["format"]:
                if format_name in sample:
                    format_column_data.append(clean_string(sample[format_name]))
                else:
                    format_column_data.append(".")

            row.extend(format_column_data)

            # headers_dict -> vcf_info_format
            info_string = make_info_string(variant)
            format_string = make_format_string(sample)

            row.append(info_string)
            row.append(format_string)

            excel_table.append(row)

    return excel_table


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

    annotated_file = vcfanno_call(vcf)
    excel_headers = make_headers()
    excel_table = make_table(annotated_file, min_vaf)

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    excel_pandasdf = pandas.DataFrame(excel_table, columns=excel_headers)

    xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
    with pandas.ExcelWriter(xlsx_name, engine="xlsxwriter") as writer:
        excel_pandasdf.to_excel(writer, sheet_name="Variants", index=False)


if __name__ == "__main__":
    do_report(
        debug=False,
        vcf="../tests/ashkenazim-trio/output3/new_normalise.vcf.gz",
        min_vaf=0.01,
        out_folder_path="../tests/ashkenazim-trio/output3",
        prefix="ashkenazim",
    )
