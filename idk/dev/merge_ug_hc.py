#!/usr/bin/env python3

"""
Based on https://git.lumc.nl/rig-framework/magpie/blob/90dd6c851959f66db944ae9915b12cc15a00de9c/scripts/merge_mity_hc.py

Merge VCF files called by GATK's HaplotypeCaller & UnifiedGenotyper


Requirements:
    * Python == 2.7.x
    * PyVCF >= 0.6.4

Copyright (c) 2013 Wibowo Arindrarto <w.arindrarto@lumc.nl>
MIT License <http://opensource.org/licenses/MIT>

"""

RELEASE = False
__version_info__ = ('0', '1', )
__version__ = '.'.join(__version_info__)
__version__ += '-dev' if not RELEASE else ''


import argparse
import os
import sys

import vcf
from vcf.parser import _Info
from vcf.utils import walk_together


def merge_hc_mity(fhc, fmity, fout, priority):
    """Merges the given HaplotypeCaller and UnifiedGenotyper VCFs into a new
    VCF."""

    hc = vcf.Reader(fhc)
    mity = vcf.Reader(fmity)

    # some sanity checks
    # TODO: possible to make it handle different samples in the two VCFs?
    if sorted(hc.samples) != sorted(mity.samples):
        raise ValueError("Input VCF files must have the same sample column headers.")
    if sorted(hc.contigs.keys()) != sorted(mity.contigs.keys()):
        raise ValueError("Input VCF files must denote the same contigs.")
    if sorted(hc.formats.keys()) != sorted(mity.formats.keys()):
        raise ValueError("Input VCF files must contain the same formats.")

    # NOTE: arbitrarily picking mity as the base template ~ we're doing
    # dict updates, so the hc values will take precedence
    # merge infos
    mity.infos.update(hc.infos)
    # merge formats ~ not necessary since they're equal
    # TODO: merge filters?
    # merge metadata
    if 'GATKCommandLine' in mity.metadata:
        mity.metadata['UnifiedGenotyperCommandLine'] = \
            mity.metadata['GATKCommandLine']
    if 'GATKCommandLine' in hc.metadata:
        mity.metadata['HaplotypeCallerCommandLine'] = \
            hc.metadata['GATKCommandLine']
    del mity.metadata['GATKCommandLine']
    del hc.metadata['GATKCommandLine']
    mity.metadata.update(hc.metadata)
    # add custom INFO field, denoting the variant caller for each variant
    # iterate over both, picking the priority when variants are called by both
    # files
    mity.infos['GATKCaller'] = _Info('GATKCaller', '.', 'String', 'GATK '
            'variant caller used to call the variant')

    out_writer = vcf.Writer(fout, mity)
    for hc_rec, mity_rec in walk_together(hc, mity):
        if hc_rec.CHROM != "MT":
          out_writer.write_record(hc_rec)
        elif mity_rec.CHROM == "MT":
          out_writer.write_record(mity_rec)
        else:
            assert False, "We should not be here!"


if __name__ == '__main__':

    usage = __doc__.split('\n\n\n')
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=usage[0], epilog=usage[1])

    def _file_read(fname):
        """Returns an open file handle if the given filename exists."""
        if not os.path.exists(fname):
            parser.error("File '{0}' not found.".format(fname))
        return open(fname, 'r')

    parser.add_argument('priority', type=str, choices=['mity', 'hc', 'both'],
            help='Which caller\'s variant is kept when a variant is called in '
            'both files')
    parser.add_argument('hc_vcf', type=_file_read, help='VCF file produced by '
            'HaplotypeCaller')
    parser.add_argument('mity_vcf', type=_file_read, help='VCF file produced by '
            'UnifiedGenotyper')
    parser.add_argument('--version', action='version', version='%(prog)s ' +
            __version__)

    args = parser.parse_args()

    merge_hc_mity(args.hc_vcf, args.mity_vcf, sys.stdout, args.priority)


hc_vcf_f = "15F00004.hc.vqsr.vcf"
mity_vcf_f = "15F00004.mity.vcf"
hc_vcf = open(hc_vcf_f, "r")
mity_vcf = open(mity_vcf_f, "r")
hc = vcf.Reader(hc_vcf)
mity = vcf.Reader(mity_vcf)
hc.contigs.keys()
mity.contigs.keys()
