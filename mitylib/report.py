import sys
import pysam
import logging
import gzip
import pandas
import os.path
import xlsxwriter
import subprocess
import yaml

# LOGGING
logger = logging.getLogger(__name__)

# def make_excel_row(variant, sample)


class Vep:
    """
    Handles everything VEP related.
    """

    def __init__(self, vcf_header) -> None:
        self.vcf_header = vcf_header
        self.vep_keys = self.get_vep_keys()

    def get_vep_dict(self, variant):
        """
        Takes two strings from VEP consequences/impacts in the form:
            impact name  | impact name  | ...
            impact value | impact value | ...

        And returns a dictionary matching relevant values.
        """

        values = variant.info["CSQ"].split("|")
        res = dict(map(lambda i, j: (i, j), self.vep_keys, values))

        return res

    def is_vepped(self, vcf):
        """
        Checks whether a file has been vepped by looking for VEP in the vcf
        header. The line should look something like this:
            ##VEP="v110" time="2023-10-04 00:20:30" cache="/net/isilonP...
        """
        if "VEP" in str(self.vcf.header):
            return True
        return False

    def get_vep_keys(self):
        """
        Example header INFO line:
            ##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|BIOTYPE|EXON|INTRON|HGVSc|HGVSp|cDNA_position|CDS_position|Protein_position|Amino_acids|Codons|Existing_variation|DISTANCE|STRAND|FLAGS|SYMBOL_SOURCE|HGNC_ID|MANE_SELECT|MANE_PLUS_CLINICAL|TSL|APPRIS|SIFT|PolyPhen|AF|CLIN_SIG|SOMATIC|PHENO|PUBMED|MOTIF_NAME|MOTIF_POS|HIGH_INF_POS|MOTIF_SCORE_CHANGE|TRANSCRIPTION_FACTORS">

        Returns a list of keys, i.e. ["Allele", "Consequence", ...]

        NOTE: This method is all hard coded, and will not work if the description format changes.
        """
        for header in self.vcf_header.records:
            if header.type == "INFO" and "CSQ" in str(header) and header.attrs:
                description = header.attrs[-2][1].strip('"')

                # change this line if the description text or format changes
                description = description.replace(
                    "Consequence annotations from Ensembl VEP. Format:", ""
                )

                keys = description.split("|")
                return keys


class Report:
    """
    TODO: Mity report description.
    """

    def __init__(self, vcf_path, min_vaf) -> None:
        self.min_vaf = min_vaf

        self.vcf_path = vcf_path
        self.vcf_obj = pysam.VariantFile(vcf_path)

        self.annot_vcf_path = None
        self.annot_vcf_obj = None

        self.vep = Vep(self.vcf_obj)

        # excel and vcf headers
        with open("report-config.yaml", "r") as report_config_file:
            self.report_config = yaml.safe_load(report_config_file)

        self.excel_headers = self.report_config["excel_headers"]
        self.vep_excel_headers = self.report_config["vep_excel_headers"]
        self.vcf_headers = self.report_config["vcf_headers"]

    def do_report(self):
        pass

    def make_hgvs(self, pos, ref, alt):
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

    def vcfanno_call(self):
        """
        Calls vcfanno to annotate the output of mity normalise with the relevant
        annotations based on vcfanno-config.toml in mitylib/
        """
        logger.debug("Running vcfanno...")

        # annotated_file name
        annotated_file = self.vcf_path.replace(".vcf.gz", ".mity.annotated.vcf")

        # vcfanno call
        vcfanno_cmd = (
            f"vcfanno -p 4 vcfanno-config.toml {self.vcf_path} > {annotated_file}"
        )
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

        self.annot_vcf_path = annotated_file
        self.annot_vcf_obj = pysam.VariantFile(annotated_file)

    def make_info_string(self, variant):
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

    def make_format_string(self, sample):
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

    def clean_string(self, s):
        """
        Removes the following characters from a string:
            "
            '
            ()

        Used to output text from annotation sources.
        """
        if isinstance(s, tuple):
            s = s[0]

        s = str(s)
        s = s.replace('"', "").replace("(", "").replace(")", "")

        return s

    def make_table(self):
        """
        Takes a vcfanno annotated vcf file and returns a formatted table with
        relevant information. The vcf names are hardcoded in vcf_column_names
        """

        excel_table = {}
        for header in self.excel_headers:
            excel_table[header] = []

        if self.vep.is_vepped():
            for header in self.vep_excel_headers:
                excel_table[header] = []

        num_samples = len(self.annot_vcf_obj.header.samples)

        for variant in self.annot_vcf_obj.fetch():
            cohort_count = 0
            info_string = self.make_info_string(variant)

            # samples
            for sample in variant.samples.values():
                # skip sample if the VAF is too low
                if float(sample["VAF"][0]) <= float(self.min_vaf):
                    continue
                
                cohort_count += 1

                excel_table["SAMPLE"] = sample.name
                excel_table["HGVS"] = self.make_hgvs(variant.pos, variant.ref, variant.alts[0])

                excel_table["CHR"] = variant.chrom
                excel_table["POS"] = variant.pos
                excel_table["REF"] = variant.ref
                excel_table["ALT"] = variant.alts[0]
                excel_table["QUAL"] = variant.qual
                excel_table["FILTER"] = variant.filter.keys()[0]

                excel_table["INFO"] = info_string
                excel_table["FORMAT"] = self.make_format_string(sample)

                for name in self.vcf_headers["info"]

                for name in vcf_column_names["start_info"]:
                    if name in variant.info.keys():
                        row.append(self.clean_string(variant.info[name]))
                    else:
                        row.append(".")

                for name in vcf_column_names["start_format"]:
                    if name in sample.keys():
                        row.append(self.clean_string(sample[name]))
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
                    info_column_data.append(self.clean_string(variant.info[field]))

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
                        row.append(self.clean_string(variant.info[name]))
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
                        row.append(self.clean_string(variant.info[name]))
                    else:
                        row.append(".")

                # headers_dict -> format
                format_column_data = []
                for format_name in vcf_column_names["format"]:
                    if format_name in sample:
                        format_column_data.append(
                            self.clean_string(sample[format_name])
                        )
                    else:
                        format_column_data.append(".")

                row.extend(format_column_data)

                # headers_dict -> vcf_info_format

            # fill in cohort count
            for _ in range(cohort_count):
                excel_table["COHORT COUNT"] = cohort_count
                excel_table["COHORT FREQUENCY"] = float(cohort_count / num_samples)

        return excel_table


def do_report(debug, vcfs, prefix=None, min_vaf=0.0, out_folder_path="."):
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

    if not os.path.exists(out_folder_path):
        os.makedirs(out_folder_path)

    xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
    with pandas.ExcelWriter(xlsx_name, engine="xlsxwriter") as writer:
        for vcf in vcfs:
            report = Report(vcf, min_vaf)
            df = report.get_df()

            sheet_name = vcf.replace("vcf.gz", "")
            df.to_excel(writer, sheet_name=sheet_name, index=False)


if __name__ == "__main__":
    do_report(
        debug=False,
        vcfs="../tests/ashkenazim-trio/output3/new_normalise.vcf.gz",
        min_vaf=0.01,
        out_folder_path="../tests/ashkenazim-trio/output4",
        prefix="ashkenazim",
    )

"""
df.insert(position, col name, list)
pandas.Dataframe(dict)
"""
