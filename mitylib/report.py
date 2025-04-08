"""
Adds annotations and generates mity excel report from VCF file.
"""

import logging
import os.path
import subprocess
from typing import Dict, Optional
import pysam
import pandas
import yaml

from vcf2pandas import vcf2pandas

from mitylib.util import MityUtil

# LOGGING
logger = logging.getLogger(__name__)


class Vep:
    """
    Provides methods to handle vepped files.
    """

    def __init__(self, vcf_obj):
        self.vcf_obj = vcf_obj
        self.vepped = self.check_vepped()

        if self.vepped:
            self.vep_keys = self.generate_vep_keys()
            self.vep_excel_headers = self.generate_vep_excel_headers()

    def get_vep_name(self, name):
        """
        Capitalises vep name, replaces _ with <space> and adds "VEP" to the end.
        """
        name = name.upper()
        name = name.replace("_", " ")
        name = name + " VEP"

        return name

    def generate_vep_excel_headers(self):
        """
        Generates vep excel headers based on the vep keys.
        """
        vep_excel_headers = ["HIGHEST IMPACT VEP"]
        for name in self.vep_keys:
            vep_excel_headers.append(self.get_vep_name(name))

        return vep_excel_headers

    def find_highest_impact(self, impacts):
        """
        Returns the highest impact in the list of consequences. Raises error if
        none of the consequences match.
        """
        if "HIGH" in impacts:
            return "HIGH"
        if "MODERATE" in impacts:
            return "MODERATE"
        if "MODIFIER" in impacts:
            return "MODIFIER"
        if "LOW" in impacts:
            return "LOW"

        logging.warning("Unknown SO term.")

        return ""

    def get_vep_values(self, variant):
        """
        Takes a string from VEP consequences/impacts in the form:
            impact value | impact value | ... |, |||, |||

        And peforms the following:
            - removes annotation "line" if the consequence has "stream" in it,
              e.g. "upstream"
            - concatenates remaining line fields with ";"

        Types:
            variant.info["CSQ"]: Tuple of form (a|b|c|..., a|b|c|..., a|b|c|...)

        Example:
            vep_keys = [ IMPACT, Consequence, field_1, field_2, field_3 ]

            list(variant.info["CSQ"] = [
                HIGH        | something_else    | a | b | c,
                LOW         | upstream_variant  | d | e | f,
                MODIFIER    | something_else    | g | h | i
            ]

            vep_dict = {
                HIGHEST IMPACT VEP  : HIGH,
                CONSEQUENCE VEP     : something_else;something_else,
                FIELD 1 VEP         : a;g,
                FIELD 2 VEP         : b;h,
                FIELD 3 VEP         : c;i
            }
        """

        vep_dict = {}
        for name in self.vep_excel_headers:
            vep_dict[name] = []

        impacts = []

        for line in list(variant.info["CSQ"]):
            split_line = line.split("|")

            # NOTE: rename IMPACT if the name changes to something else
            impact = split_line[self.vep_keys.index("IMPACT")]
            if "stream" in impact:
                continue

            impacts.append(impact)

            for i, value in enumerate(split_line):
                if value:
                    vep_dict[self.vep_excel_headers[i + 1]].append(value)

        vep_dict["HIGHEST IMPACT VEP"] = self.find_highest_impact(impacts)

        # NOTE: final string format/join for excel output can be changed here
        for key, value in vep_dict.items():
            if isinstance(value, list):
                vep_dict[key] = ";".join(value)

        return vep_dict

    def check_vepped(self):
        """
        Checks whether a file has been vepped by looking for VEP in the vcf
        header. The line should look something like this:
            ##VEP="v110" time="2023-10-04 00:20:30" cache="/net/isilonP...
        """
        if "VEP" in str(self.vcf_obj.header):
            return True
        return False

    def generate_vep_keys(self):
        """
        Example header INFO line:
            ##INFO=<ID=CSQ,Number=.,Type=String,Description=
            "Consequence annotations from Ensembl VEP. Format:
            Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...

        Returns a list of keys, i.e. ["Allele", "Consequence", ...]

        NOTE: This method is hard coded based on the format of the description.
        """
        for header in self.vcf_obj.header.records:
            if header.type == "INFO" and "CSQ" in str(header) and header.attrs:
                description = header.attrs[-2][1].strip('"')

                # change this line if the description text or format changes
                description = description.replace(
                    "Consequence annotations from Ensembl VEP. Format:", ""
                )

                keys = description.split("|")
                return keys

        return None


class SingleReport:
    """
    Handles generating a mity report for one VCF file.
    """

    def __init__(
        self,
        vcf_path: str,
        min_vaf: str,
        keep: bool,
        vcfanno_base_path: str,
        vcfanno_config: str,
        report_config: str,
        prefix: str,
        output_dir: str,
        output_annotated_vcf: bool = False,
    ) -> None:
        self.min_vaf = min_vaf
        self.keep = keep

        self.vcf_path = vcf_path
        self.vcf_obj = pysam.VariantFile(vcf_path)

        self.annot_vcf_path: Optional[str] = None
        self.annot_vcf_obj = None

        self.vcfanno_base_path = vcfanno_base_path
        self.vcfanno_config = vcfanno_config
        self.report_config = report_config

        self.prefix = prefix
        self.output_dir = output_dir
        self.output_annotated_vcf = output_annotated_vcf

        # excel and vcf headers
        with open(self.report_config, "r", encoding="utf-8") as report_config_file:
            opened_report_config = yaml.safe_load(report_config_file)

        self.excel_headers = opened_report_config["excel_headers"]
        self.vcf_headers = opened_report_config["vcf_headers"]

        self.excel_table: Dict[str, list[str]] = {}
        self.df = None

        self.vep: Optional[Vep] = None

        self.run()

    def run(self):
        """
        Run SingleReport.
        """

        self.vcfanno_call()
        self.vep = Vep(self.annot_vcf_obj)
        self.make_table()

        if self.output_annotated_vcf:
            annotated_vcf_df = vcf2pandas(self.annot_vcf_path)
            annotated_vcf_df.to_excel(
                os.path.join(self.output_dir, self.prefix + ".mity.annotated.vcf.xlsx"),
                index=False,
            )

        if not self.keep:
            # remove vcfanno annotated vcf
            os.remove(self.annot_vcf_path)

    def get_df(self):
        """
        Return a pandas dataframe from the generated excel table.
        """
        all_excel_headers = self.excel_headers
        if self.vep.vepped:
            all_excel_headers += self.vep.vep_excel_headers
        return pandas.DataFrame(self.excel_table, columns=all_excel_headers)

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
        base_path_arg = f"-base-path {self.vcfanno_base_path}" if self.vcfanno_base_path else ""

        # vcfanno call
        vcfanno_cmd = (
            f"vcfanno -p 4 {base_path_arg} {self.vcfanno_config} {self.vcf_path} > {annotated_file}"
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

    def clean_string(self, s: str) -> str:
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
        Takes a vcfanno annotated vcf file and returns a formatted dictionary
        with relevant information.

        The vcf and excel header names are hardcoded in report-config.yaml.
        """

        for header in self.excel_headers:
            self.excel_table[header] = []

        # check if VCF file is vepped, adds relevant VEP headers

        if self.vep.vepped:
            for header in self.vep.vep_excel_headers:
                self.excel_table[header] = []

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

                self.excel_table["SAMPLE"].append(sample.name)
                self.excel_table["HGVS"].append(
                    self.make_hgvs(variant.pos, variant.ref, variant.alts[0])
                )

                self.excel_table["CHR"].append(variant.chrom)
                self.excel_table["POS"].append(variant.pos)
                self.excel_table["REF"].append(variant.ref)

                # NOTE: assumes only one ALT as the vcf should be normalised
                self.excel_table["ALT"].append(variant.alts[0])
                self.excel_table["QUAL"].append(variant.qual)
                self.excel_table["FILTER"].append(variant.filter.keys()[0])

                self.excel_table["INFO"].append(info_string)
                self.excel_table["FORMAT"].append(self.make_format_string(sample))

                # vcf_headers: info
                for vcf_header, excel_header in self.vcf_headers["info"].items():
                    if vcf_header in variant.info.keys():
                        self.excel_table[excel_header].append(
                            self.clean_string(variant.info[vcf_header])
                        )
                    else:
                        self.excel_table[excel_header].append(".")

                # vcf_headers: annotations
                for vcf_header, excel_header in self.vcf_headers["annotations"].items():
                    if vcf_header in variant.info.keys():
                        self.excel_table[excel_header].append(
                            self.clean_string(variant.info[vcf_header])
                        )
                    else:
                        self.excel_table[excel_header].append(".")

                # vcf_headers: format
                for vcf_header, excel_header in self.vcf_headers["format"].items():
                    if vcf_header in sample.keys():
                        self.excel_table[excel_header].append(self.clean_string(sample[vcf_header]))
                    else:
                        self.excel_table[excel_header].append(".")

            # fill in cohort count
            cohort_frequency = float(cohort_count / num_samples)
            for _ in range(cohort_count):
                self.excel_table["COHORT COUNT"].append(cohort_count)
                self.excel_table["COHORT FREQUENCY"].append(cohort_frequency)

            # fill in VEP impacts if vepped
            if self.vep.vepped:
                self.add_vep_impacts(variant, cohort_count)

    def add_vep_impacts(self, variant, cohort_count: int) -> None:
        """
        Adds vep impacts for a variant.
        """
        if not self.vep:
            raise RuntimeError("Vep not initialised, something likely went wrong with vcfanno")

        variant_impacts = self.vep.get_vep_values(variant)
        for _ in range(cohort_count):
            for vep_header in self.vep.vep_excel_headers:
                self.excel_table[vep_header].append(variant_impacts[vep_header])


class Report:
    """
    Runs mity report.
    """

    def __init__(
        self,
        debug: bool,
        vcfs,
        contig: str,
        prefix: Optional[str] = None,
        min_vaf: float = 0.0,
        output_dir: str = ".",
        keep: bool = False,
        vcfanno_base_path: Optional[str] = None,
        vcfanno_config: Optional[str] = None,
        report_config: Optional[str] = None,
        output_annotated_vcf: bool = False,
    ) -> None:
        self.debug = debug
        self.vcfs = vcfs[0]
        self.contig = contig
        self.prefix = prefix
        self.min_vaf = min_vaf
        self.output_dir = output_dir
        self.keep = keep
        self.vcfanno_base_path = vcfanno_base_path
        self.vcfanno_config = vcfanno_config
        self.report_config = report_config
        self.output_annotated_vcf = output_annotated_vcf

        self.run()

    def run(self) -> None:
        """
        Runs mity report.
        """
        if self.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Entered debug mode.")
        else:
            logger.setLevel(logging.INFO)

        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        # self.vcfs is of type str if there is only one vcf input, but list[str]
        # so we standardise everything to be a list[str]
        if isinstance(self.vcfs, str):
            self.vcfs = [self.vcfs]

        if self.prefix is None:
            self.prefix = MityUtil.make_prefix(self.vcfs[0])

        config_path = os.path.join(MityUtil.get_mity_dir(), "config")
        if self.vcfanno_config is None:
            self.vcfanno_base_path = MityUtil.get_mity_dir()
            match self.contig:
                case "MT":
                    self.vcfanno_config = os.path.join(config_path, "vcfanno-config-mt.toml")
                case "chrM":
                    self.vcfanno_config = os.path.join(config_path, "vcfanno-config-chrm.toml")
                case _:
                    raise ValueError(
                        "Contig not recognised, please specify a valid contig (either MT or chrM)"
                    )

        if self.report_config is None:
            self.report_config = os.path.join(config_path, "report-config.yaml")

        xlsx_name = os.path.join(self.output_dir, self.prefix + ".mity.report.xlsx")
        with pandas.ExcelWriter(xlsx_name, engine="xlsxwriter") as writer:
            for vcf in self.vcfs:
                single_report = SingleReport(
                    vcf_path=vcf,
                    min_vaf=self.min_vaf,
                    keep=self.keep,
                    vcfanno_base_path=self.vcfanno_base_path,
                    vcfanno_config=self.vcfanno_config,
                    report_config=self.report_config,
                    prefix=self.prefix,
                    output_dir=self.output_dir,
                    output_annotated_vcf=self.output_annotated_vcf,
                )
                df = single_report.get_df()

                sheet_name = vcf.replace(".vcf.gz", "").split("/")[-1]
                if len(sheet_name) > 31:
                    logging.info(
                        "sheet_name: %s was too long and was automatically shortened",
                        sheet_name,
                    )
                sheet_name = sheet_name[:31]
                df.to_excel(writer, sheet_name=sheet_name, index=False)
