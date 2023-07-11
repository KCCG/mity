"""
Mity: a sensitive variant analysis pipeline optimised for WGS data



Usage: See the online manual for details: http://github.com/KCCG/mity
Authors: Clare Puttick, Mark Cowley
License: MIT
"""
import argparse
import logging

from . import (call, normalise, report, merge)
from ._version import __version__
from .util import select_reference_fasta
from .util import select_reference_genome

__all__ = []


def public(fn):
    __all__.append(fn.__name__)
    return fn


usage = __doc__.split('\n\n\n', maxsplit=1)
usage[-1] += "Version: " + __version__

AP = argparse.ArgumentParser(description=usage[0], epilog=usage[1], formatter_class=argparse.RawTextHelpFormatter)
AP_subparsers = AP.add_subparsers(
        help="mity sub-commands (use with -h for more info)")

# call -------------------------------------------------------------------------

do_call = public(call.do_call)

def _cmd_call(args):
    if True:
        # TODO: why does this not turn on debugging output? It needs to be turned on in `mity:L20` to work
        logging.basicConfig(level=logging.DEBUG, format="%(message)s")
        #logging.getLogger('mitylib').setLevel(level=logging.DEBUG)

    """Call mitochondrial variants"""
    logging.info("mity version %s", __version__)
    logging.info("Calling mitochondrial variants")
    logging.debug("Debugging mode activated")

    genome = select_reference_genome(args.reference, None)
    args.reference = select_reference_fasta(args.reference, None)

    call.do_call(args.bam, args.reference, genome, args.prefix, args.min_mq,
                 args.min_bq, args.min_af, args.min_ac, args.p, args.normalise,
                 args.out_folder_path, args.region)


P_call = AP_subparsers.add_parser('call', help=_cmd_call.__doc__)
P_call.add_argument('bam', action='append', nargs='+',
                    help='BAM files to run the analysis on.')
P_call.add_argument('--reference', choices=['hs37d5', 'hg19', 'hg38', 'mm10'],
                    default="hs37d5", required=False,
                    help='reference genome version to use. default: hs37d5')
# P_call.add_argument('--custom_reference', action='store',
#                     default="", required=False,
#                     help='The path to a custom reference genome file in uncompressed fasta format')
P_call.add_argument('--prefix', action='store',
                    help='Output files will be named with PREFIX')
P_call.add_argument('--min-mapping-quality', action='store', type=int,
                    default=30,
                    help='Exclude alignments from analysis if they have a '
                         'mapping quality less than MIN_MAPPING_QUALITY. '
                         'Default: 30',
                    dest="min_mq")
P_call.add_argument('--min-base-quality', action='store', type=int, default=24,
                    help='Exclude alleles from analysis if their supporting '
                         'base quality is less than MIN_BASE_QUALITY. '
                         'Default: 24',
                    dest="min_bq")
P_call.add_argument('--min-alternate-fraction', action='store', type=float,
                    default=0.01,
                    help='Require at least MIN_ALTERNATE_FRACTION '
                         'observations supporting an alternate allele within '
                         'a single individual in the in order to evaluate the '
                         'position. Default: 0.01, range = [0,1]',
                    dest="min_af")
P_call.add_argument('--min-alternate-count', action='store', type=int,
                    default=4,
                    help='Require at least MIN_ALTERNATE_COUNT observations '
                         'supporting an alternate allele within a single '
                         'individual in order to evaluate the position. '
                         'Default: 4',
                    dest="min_ac")
P_call.add_argument('--p', action='store', type=float,
                    default=0.002,
                    help='Minimum noise level. This is used to calculate QUAL score. '
                         'Default: 0.002, range = [0,1]',
                    dest="p")
P_call.add_argument('--normalise', action='store_true',
                    help='Normalise the resulting VCF?')
P_call.add_argument('--out-folder-path', action='store', type=str,
                    default='.',
                    help='Output files will be saved in OUT_FOLDER_PATH. '
                         "Default: '.' ",
                    dest="out_folder_path")
P_call.add_argument('--region', action='store', type=str,
                    default=None,
                    help='Region of MT genome to call variants in. '
                         'If unset will call variants in entire MT genome as specified in BAM header. '
                         "Default: Entire MT genome. ",
                    dest="region")
P_call.add_argument('--debug', action='store_true',
                    help='Verbose output for debugging?')
P_call.add_argument('--bam-file-list', action='store', type=str,
                    default=None,
                    help='A text file of BAM files to be processed. The path to each file should be on one row per Region of MT genome to call variants in. '
                         'If unset will call variants in entire MT genome as specified in BAM header. '
                         "Default: Entire MT genome. ",
                    dest="region")
P_call.set_defaults(func=_cmd_call)

# normalise --------------------------------------------------------------------

do_normalise = public(normalise.do_normalise)


def _cmd_normalise(args):
    """Normalise & FILTER mitochondrial variants"""
    logging.info("mity %s", __version__)
    logging.info("Normalising and FILTERing mitochondrial vcf.gz file")

    genome = select_reference_genome(args.reference, None)

    normalise.do_normalise(args.vcf, args.outfile, p=args.p, SB_range=[0.1, 0.9], genome=genome)


P_normalise = AP_subparsers.add_parser('normalise', help=_cmd_normalise.__doc__)
P_normalise.add_argument('vcf', action='store',
                         help="vcf.gz file from running mity")
P_normalise.add_argument('--outfile', action='store', required=True,
                         help="output VCF file in bgzip compressed format")
P_normalise.add_argument('--p', action='store', type=float,
                         default=0.002,
                         help='Minimum noise level. This is used to calculate QUAL score'
                              'Default: 0.002, range = [0,1]',
                         dest="p")
P_normalise.add_argument('--reference', choices=['hs37d5', 'hg19', 'hg38', 'mm10'],
                    default="hs37d5", required=False,
                    help='reference genome version to use. default: hs37d5')
P_normalise.set_defaults(func=_cmd_normalise)


# report -----------------------------------------------------------------------

do_report = public(report.do_report)


def _cmd_report(args):
    """Generate mity report"""
    logging.info("mity %s", __version__)
    logging.info("Generating mity report")
    report.do_report(args.vcf, args.prefix, args.min_vaf, args.out_folder_path)


P_report = AP_subparsers.add_parser('report', help=_cmd_report.__doc__)
P_report.add_argument('--prefix', action='store',
                      help='Output files will be named with PREFIX')
P_report.add_argument('--min_vaf', action='store', type=float, default=0, help=
'A variant must have at least this VAF to be included in the report. Default: '
'0.')
P_report.add_argument('--out-folder-path', action='store', type=str,
                    default='.',
                    help='Output files will be saved in OUT_FOLDER_PATH. '
                         "Default: '.' ",
                    dest="out_folder_path")
P_report.add_argument('vcf', action='append', nargs='+',
                    help="mity vcf files to create a report from")
P_report.set_defaults(func=_cmd_report)


# merge -----------------------------------------------------------------------

do_merge = public(merge.do_merge)

def _cmd_merge(args):
    """Merging mity VCF with nuclear VCF"""
    logging.info("mity %s", __version__)
    logging.info("mity vcf merge")

    genome = select_reference_genome(args.reference, None)

    merge.do_merge(args.mity_vcf, args.nuclear_vcf, args.prefix, genome)

P_merge = AP_subparsers.add_parser('merge', help=_cmd_merge.__doc__)
P_merge.add_argument('--mity_vcf', action='store', required=True,
                      help="mity vcf file")
P_merge.add_argument('--nuclear_vcf', action='store', required=True,
                     help="nuclear vcf file")
P_merge.add_argument('--prefix', action='store',
                     help='Output files will be named with PREFIX. '
                     'The default is to use the nuclear vcf name')
P_merge.add_argument('--reference', choices=['hs37d5', 'hg19', 'hg38', 'mm10'],
                     default="hs37d5", required=False,
                     help='reference genome version to use. default: hs37d5')
# P_merge.add_argument('--custom_reference', action='store',
#                      default="", required=False,
#                      help='The path to a custom reference genome file in uncompressed fasta format')
P_merge.set_defaults(func=_cmd_merge)


# version ----------------------------------------------------------------------

def print_version(_args):
    """Display this program's version."""
    print(__version__)


P_version = AP_subparsers.add_parser('version', help=print_version.__doc__)
P_version.set_defaults(func=print_version)


def parse_args(args=None):
    """Parse the command line."""
    return AP.parse_args(args=args)
