"""
Mity: a sensitive variant analysis pipeline optimised for WGS data



Usage: See the online manual for details: http://github.com/KCCG/mity
Authors: Clare Puttick, Mark Cowley, Trent Zeng, Christian Fares
License: MIT
"""

import argparse
import logging
import os

from mitylib import call, normalise, report, merge, util
from ._version import __version__


usage = __doc__.split("\n\n\n", maxsplit=1)
usage[-1] += "Version: " + __version__

AP = argparse.ArgumentParser(
    description=usage[0], epilog=usage[1], formatter_class=argparse.RawTextHelpFormatter
)

# subparser    -----------------------------------------------------------------

AP_subparsers = AP.add_subparsers(help="mity sub-commands (use with -h for more info)")

# call -------------------------------------------------------------------------


def _cmd_call(args):
    """Call mitochondrial variants"""
    logging.info("mity version %s", __version__)
    logging.info("Calling mitochondrial variants")
    logging.debug("Debugging mode activated")

    genome = util.MityUtil.select_reference_genome(args.reference, args.custom_reference_genome)
    args.reference = util.MityUtil.select_reference_fasta(args.reference, args.custom_reference_fasta)

    call.Call(
        debug=args.debug,
        files=args.files,
        reference=args.reference,
        genome=genome,
        prefix=args.prefix,
        min_mq=args.min_mq,
        min_bq=args.min_bq,
        min_af=args.min_af,
        min_ac=args.min_ac,
        p=args.p,
        normalise=args.normalise,
        output_dir=args.output_dir,
        region=args.region,
        bam_list=args.bam_file_list,
        keep=args.keep,
    )


P_call = AP_subparsers.add_parser("call", help=_cmd_call.__doc__)
P_call.add_argument(
    "-d", "--debug", action="store_true", help="Enter debug mode", required=False
)
P_call.add_argument(
    "files",
    action="append",
    nargs="+",
    help="BAM / CRAM files to run the analysis on. If --bam-file-list is included, this argument is the file containing the list of bam/cram files.",
)
P_call.add_argument(
    "--reference",
    choices=["hs37d5", "hg19", "hg38", "mm10"],
    default="hs37d5",
    required=False,
    help="Reference genome version to use. Default: hs37d5",
)
P_call.add_argument(
    "--custom-reference-fasta",
    action="store",
    help="Specify custom reference fasta file",
    dest="custom_reference_fasta"
)
P_call.add_argument(
    "--custom-reference-genome",
    action="store",
    help="Specify custom reference genome file",
    dest="custom_reference_genome"
)
P_call.add_argument(
    "--prefix", action="store", help="Output files will be named with PREFIX"
)
P_call.add_argument(
    "--min-mapping-quality",
    action="store",
    type=int,
    default=30,
    help="Exclude alignments from analysis if they have a "
    "mapping quality less than MIN_MAPPING_QUALITY. "
    "Default: 30",
    dest="min_mq",
)
P_call.add_argument(
    "--min-base-quality",
    action="store",
    type=int,
    default=24,
    help="Exclude alleles from analysis if their supporting "
    "base quality is less than MIN_BASE_QUALITY. "
    "Default: 24",
    dest="min_bq",
)
P_call.add_argument(
    "--min-alternate-fraction",
    action="store",
    type=float,
    default=0.01,
    help="Require at least MIN_ALTERNATE_FRACTION "
    "observations supporting an alternate allele within "
    "a single individual in the in order to evaluate the "
    "position. Default: 0.01, range = [0,1]",
    dest="min_af",
)
P_call.add_argument(
    "--min-alternate-count",
    action="store",
    type=int,
    default=4,
    help="Require at least MIN_ALTERNATE_COUNT observations "
    "supporting an alternate allele within a single "
    "individual in order to evaluate the position. "
    "Default: 4",
    dest="min_ac",
)
P_call.add_argument(
    "--p",
    action="store",
    type=float,
    default=0.002,
    help="Minimum noise level. This is used to calculate QUAL score. "
    "Default: 0.002, range = [0,1]",
    dest="p",
)
P_call.add_argument(
    "--normalise", action="store_true", help="Run mity normalise the resulting VCF"
)
P_call.add_argument(
    "--output-dir",
    action="store",
    type=str,
    default=".",
    help="Output files will be saved in OUTPUT_DIR. " "Default: '.' ",
    dest="output_dir",
)
P_call.add_argument(
    "--region",
    action="store",
    type=str,
    default=None,
    help="Region of MT genome to call variants in. "
    "If unset will call variants in entire MT genome as specified in BAM header. "
    "Default: Entire MT genome. ",
    dest="region",
)
P_call.add_argument(
    "--bam-file-list",
    action="store_true",
    default=False,
    help="Treat the file as a text file of BAM files to be processed."
    " The path to each file should be on one row per bam file.",
    dest="bam_file_list",
)
P_call.add_argument(
    "-k",
    "--keep",
    action="store_true",
    required=False,
    help="Keep all intermediate files",
)
P_call.set_defaults(func=_cmd_call)

# normalise --------------------------------------------------------------------


def _cmd_normalise(args):
    """Normalise & FILTER mitochondrial variants"""
    logging.info("mity %s", __version__)
    logging.info("Normalising and FILTERing mitochondrial vcf.gz file")

    genome = util.MityUtil.select_reference_genome(args.reference, args.custom_reference_genome)
    args.reference = util.MityUtil.select_reference_fasta(args.reference, args.custom_reference_fasta)

    normalise.Normalise(
        debug=args.debug,
        vcf=args.vcf,
        reference_fasta=args.reference,
        genome=genome,
        output_dir=args.output_dir,
        prefix=args.prefix,
        allsamples=args.allsamples,
        keep=args.keep,
        p=args.p,
    )


P_normalise = AP_subparsers.add_parser("normalise", help=_cmd_normalise.__doc__)
P_normalise.add_argument(
    "-d", "--debug", action="store_true", help="Enter debug mode", required=False
)
P_normalise.add_argument("vcf", action="store", help="vcf.gz file from running mity")
P_normalise.add_argument(
    "--output-dir",
    action="store",
    type=str,
    default=".",
    help="Output files will be saved in OUTPUT_DIR. " "Default: '.' ",
    dest="output_dir",
)
P_normalise.add_argument(
    "--prefix", action="store", help="Output files will be named with PREFIX"
)
P_normalise.add_argument(
    "--allsamples",
    action="store_true",
    required=False,
    help="PASS in the filter requires all samples to pass instead of just one",
)
P_normalise.add_argument(
    "-k",
    "--keep",
    action="store_true",
    required=False,
    help="Keep all intermediate files",
)
P_normalise.add_argument(
    "--p",
    action="store",
    type=float,
    default=0.002,
    help="Minimum noise level. This is used to calculate QUAL score"
    "Default: 0.002, range = [0,1]",
    dest="p",
)
P_normalise.add_argument(
    "--reference",
    choices=["hs37d5", "hg19", "hg38", "mm10"],
    default="hs37d5",
    required=False,
    help="Reference genome version to use. default: hs37d5",
)
P_normalise.add_argument(
    "--custom-reference-fasta",
    action="store",
    help="Specify custom reference fasta file",
    dest="custom_reference_fasta"
)
P_normalise.add_argument(
    "--custom-reference-genome",
    action="store",
    help="Specify custom reference genome file",
    dest="custom_reference_genome"
)
P_normalise.set_defaults(func=_cmd_normalise)


# report -----------------------------------------------------------------------


def _cmd_report(args):
    """Generate mity report"""
    logging.info("mity %s", __version__)
    logging.info("Generating mity report")
    report.Report(
        debug=args.debug,
        vcfs=args.vcf,
        contig=args.contig,
        prefix=args.prefix,
        min_vaf=args.min_vaf,
        output_dir=args.output_dir,
        keep=args.keep,
        vcfanno_base_path=args.vcfanno_base_path,
        vcfanno_config=args.vcfanno_config,
        report_config=args.report_config,
        output_annotated_vcf=args.output_annotated_vcf,
    )


P_report = AP_subparsers.add_parser("report", help=_cmd_report.__doc__)
P_report.add_argument(
    "-d", "--debug", action="store_true", help="Enter debug mode", required=False
)
P_report.add_argument(
    "--prefix", action="store", help="Output files will be named with PREFIX"
)
P_report.add_argument(
    "--min_vaf",
    action="store",
    type=float,
    default=0,
    help="A variant must have at least this VAF to be included in the report. Default: "
    "0.",
)
P_report.add_argument(
    "--output-dir",
    action="store",
    type=str,
    default=".",
    help="Output files will be saved in OUTPUT_DIR. " "Default: '.' ",
    dest="output_dir",
)
P_report.add_argument(
    "vcf", action="append", nargs="+", help="mity vcf files to create a report from"
)
P_report.add_argument(
    "-k",
    "--keep",
    action="store_true",
    required=False,
    help="Keep all intermediate files",
)
P_report.add_argument(
    "--contig",
    choices=["MT", "chrM"],
    default="MT",
    required=False,
    help="Contig used for annotation purposes",
)
P_report.add_argument(
    "--vcfanno-base-path",
    action="store",
    help="Path to the custom annotations used for vcfanno. Only required if using custom annotations.",
    dest="vcfanno_base_path",
)
P_report.add_argument(
    "--custom-vcfanno-config",
    action="store",
    help="Provide a custom vcfanno-config.toml for custom annotations.",
    dest="vcfanno_config",
)
P_report.add_argument(
    "--custom-report-config",
    action="store",
    help="Provide a custom report-config.yaml for custom report generation.",
    dest="report_config",
)
P_report.add_argument(
    "--output-annotated-vcf",
    action="store_true",
    help="Output annotated vcf file",
    dest="output_annotated_vcf",
)

P_report.set_defaults(func=_cmd_report)


# merge -----------------------------------------------------------------------


def _cmd_merge(args):
    """Merging mity VCF with nuclear VCF"""
    logging.info("mity %s", __version__)
    logging.info("mity vcf merge")

    genome = util.MityUtil.select_reference_genome(args.reference, None)

    merge.Merge(
        debug=args.debug,
        nuclear_vcf_path=args.nuclear_vcf,
        mity_vcf_path=args.mity_vcf,
        genome=genome,
        output_dir=args.output_dir,
        prefix=args.prefix,
        keep=args.keep,
    )


P_merge = AP_subparsers.add_parser("merge", help=_cmd_merge.__doc__)
P_merge.add_argument("--mity_vcf", action="store", required=True, help="mity vcf file")
P_merge.add_argument(
    "--nuclear_vcf", action="store", required=True, help="nuclear vcf file"
)
P_merge.add_argument(
    "--output-dir",
    action="store",
    type=str,
    default=".",
    help="Output files will be saved in OUTPUT_DIR. " "Default: '.' ",
    dest="output_dir",
)
P_merge.add_argument(
    "--prefix",
    action="store",
    help="Output files will be named with PREFIX. "
    "The default is to use the nuclear vcf name",
)
P_merge.add_argument(
    "--reference",
    choices=["hs37d5", "hg19", "hg38", "mm10"],
    default="hs37d5",
    required=False,
    help="reference genome version to use. default: hs37d5",
)
P_merge.add_argument(
    "-d", "--debug", action="store_true", help="Enter debug mode", required=False
)
P_merge.add_argument(
    "-k",
    "--keep",
    action="store_true",
    required=False,
    help="Keep all intermediate files",
)
P_merge.set_defaults(func=_cmd_merge)


# runall -----------------------------------------------------------------------


def _cmd_runall(args):
    """Run MITY call, normalise and report all in one go."""
    logging.info("mity %s", __version__)
    logging.info("mity runall")

    genome = util.MityUtil.select_reference_genome(args.reference, args.custom_reference_genome)
    args.reference = util.MityUtil.select_reference_fasta(args.reference, args.custom_reference_fasta)

    call.Call(
        debug=args.debug,
        files=args.files,
        reference=args.reference,
        genome=genome,
        prefix=args.prefix,
        min_mq=args.min_mq,
        min_bq=args.min_bq,
        min_af=args.min_af,
        min_ac=args.min_ac,
        p=args.p,
        # normalise flag is set to true instead of running normalise separately
        normalise=True,
        output_dir=args.output_dir,
        region=args.region,
        bam_list=args.bam_file_list,
        keep=args.keep,
    )

    logging.debug("mity call and normalise completed")

    # This makes use of the uniform naming scheme of mity command outputs.
    normalised_vcf_path = os.path.join(
        args.output_dir, args.prefix + ".normalise.vcf.gz"
    )

    logging.debug("assumed mity normalise vcf output path is: %s", normalised_vcf_path)

    # matching argparse quirk
    normalised_vcf_path = [normalised_vcf_path]

    report.Report(
        debug=args.debug,
        vcfs=normalised_vcf_path,
        contig=args.contig,
        prefix=args.prefix,
        min_vaf=args.min_vaf,
        output_dir=args.output_dir,
        keep=args.keep,
        vcfanno_base_path=args.vcfanno_base_path,
        vcfanno_config=args.vcfanno_config,
        report_config=args.report_config,
        output_annotated_vcf=args.output_annotated_vcf,
    )


P_runall = AP_subparsers.add_parser("runall", help=_cmd_runall.__doc__)
P_runall.add_argument(
    "-d", "--debug", action="store_true", help="Enter debug mode", required=False
)
P_runall.add_argument(
    "files",
    action="append",
    nargs="+",
    help="BAM / CRAM files to run the analysis on. If --bam-file-list is included, this argument is the file containing the list of bam/cram files.",
)
P_runall.add_argument(
    "--reference",
    choices=["hs37d5", "hg19", "hg38", "mm10"],
    default="hs37d5",
    required=False,
    help="Reference genome version to use. Default: hs37d5",
)
# For the runall command, we mandate that the prefix option is set. This is not
# true for regular mity call, normalise or report separately.
P_runall.add_argument(
    "--prefix",
    action="store",
    required=True,
    help="Output files will be named with PREFIX",
)
P_runall.add_argument(
    "--min-mapping-quality",
    action="store",
    type=int,
    default=30,
    help="Exclude alignments from analysis if they have a "
    "mapping quality less than MIN_MAPPING_QUALITY. "
    "Default: 30",
    dest="min_mq",
)
P_runall.add_argument(
    "--min-base-quality",
    action="store",
    type=int,
    default=24,
    help="Exclude alleles from analysis if their supporting "
    "base quality is less than MIN_BASE_QUALITY. "
    "Default: 24",
    dest="min_bq",
)
P_runall.add_argument(
    "--min-alternate-fraction",
    action="store",
    type=float,
    default=0.01,
    help="Require at least MIN_ALTERNATE_FRACTION "
    "observations supporting an alternate allele within "
    "a single individual in the in order to evaluate the "
    "position. Default: 0.01, range = [0,1]",
    dest="min_af",
)
P_runall.add_argument(
    "--min-alternate-count",
    action="store",
    type=int,
    default=4,
    help="Require at least MIN_ALTERNATE_COUNT observations "
    "supporting an alternate allele within a single "
    "individual in order to evaluate the position. "
    "Default: 4",
    dest="min_ac",
)
P_runall.add_argument(
    "--p",
    action="store",
    type=float,
    default=0.002,
    help="Minimum noise level. This is used to calculate QUAL score. "
    "Default: 0.002, range = [0,1]",
    dest="p",
)
P_runall.add_argument(
    "--output-dir",
    action="store",
    type=str,
    default=".",
    help="Output files will be saved in OUTPUT_DIR. " "Default: '.' ",
    dest="output_dir",
)
P_runall.add_argument(
    "--region",
    action="store",
    type=str,
    default=None,
    help="Region of MT genome to call variants in. "
    "If unset will call variants in entire MT genome as specified in BAM header. "
    "Default: Entire MT genome. ",
    dest="region",
)
P_runall.add_argument(
    "--bam-file-list",
    action="store_true",
    default=False,
    help="Treat the file as a text file of BAM files to be processed."
    " The path to each file should be on one row per bam file.",
    dest="bam_file_list",
)
P_runall.add_argument(
    "-k",
    "--keep",
    action="store_true",
    required=False,
    help="Keep all intermediate files",
)
P_runall.add_argument(
    "--min_vaf",
    action="store",
    type=float,
    default=0,
    help="A variant must have at least this VAF to be included in the report. Default: "
    "0.",
)
P_runall.add_argument(
    "--contig",
    choices=["MT", "chrM"],
    default="MT",
    required=False,
    help="Contig used for annotation purposes",
)
P_runall.add_argument(
    "--vcfanno-base-path",
    action="store",
    help="Path to the custom annotations used for vcfanno. Only required if using custom annotations.",
    dest="vcfanno_base_path",
)
P_runall.add_argument(
    "--custom-vcfanno-config",
    action="store",
    help="Provide a custom vcfanno-config.toml for custom annotations.",
    dest="vcfanno_config",
)
P_runall.add_argument(
    "--custom-report-config",
    action="store",
    help="Provide a custom report-config.yaml for custom report generation.",
    dest="report_config",
)
P_runall.add_argument(
    "--custom-reference-fasta",
    action="store",
    help="Specify custom reference fasta file",
    dest="custom_reference_fasta"
)
P_runall.add_argument(
    "--custom-reference-genome",
    action="store",
    help="Specify custom reference genome file",
    dest="custom_reference_genome"
)
P_runall.add_argument(
    "--output-annotated-vcf",
    action="store_true",
    help="Output annotated vcf file",
    dest="output_annotated_vcf",
)
P_runall.set_defaults(func=_cmd_runall)


# version ----------------------------------------------------------------------


def print_version(_args):
    """Display this program's version."""
    print(__version__)


P_version = AP_subparsers.add_parser("version", help=print_version.__doc__)
P_version.set_defaults(func=print_version)

# parse_args -------------------------------------------------------------------


def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
