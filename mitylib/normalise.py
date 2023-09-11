"""Normalise mity VCF."""

import subprocess
import logging
from .util import tabix
from .util import get_mity_dir
from scipy.stats import binom
from math import isinf
import numpy as np
import configparser
import pysam
import pysam.bcftools
import os.path

# CONSTANTS
P_VAL = 0.002
SB_RANGE_LO = 0.1 
SB_RANGE_HI = 0.9
MIN_MQMR = 30
MIN_AQR = 20
MIN_DP = 15

BLACKLIST = list(range(302, 319)) + list(range(3105,3109))

# CONFIG PARSER
config = configparser.ConfigParser()
config.read(get_mity_dir() + "/config.ini")
GENOME_FILE = config.get('PATHS', 'GENOME_FILE')

# LOGGING
logger = logging.getLogger(__name__)


def add_headers(input_vcf):
    # setting headers
    new_header = input_vcf.header

    # filter headers
    new_header.filters.add("POS", None, None, "Variant falls in the blacklist of positions: MT:302-319, MT:3105-3108")
    new_header.filters.add("SBR", None, None, "For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)")
    new_header.filters.add("SBA", None, None, "For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)")
    new_header.filters.add("MQMR", None, None, "For all alleles MQMR<30")
    new_header.filters.add("AQR", None, None, "For all alleles AQR<20")

    # info headers
    new_header.info.add("SBR", number=1, type="Float", description="For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)")
    new_header.info.add("SBA", number=1, type="Float", description="For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)")
    
    # format headers
    new_header.formats.add("AQA", number="A", type="Float", description="Average base quality of the alternate reads, AQA=QA/AO")
    new_header.formats.add("AQR", number="A", type="Float", description="Average base quality of the reference reads, AQR=QR/RO")
    new_header.formats.add("VAF", number="A", type="Float", description="Allele frequency in the range (0,1] - the ratio of the number of alternate reads to reference reads")
    new_header.formats.add("q", number="A", type="Float", description="Phred scaled binomial probability of seeing AO reads from DP, assuming a noise floor of p=0.002. ")
    new_header.formats.add("tier", number="A", type="Integer", description="Custom variant tier (TODO: change description)")


    return new_header

def mity_qual(AO, DP, p):
    """
    Compute variant quality
    :param AO: (int) number of alternative reads
    :param DP: (int) total read depth
    :return: (float) phred-scaled quality score

    >>> [mity_qual(x,10) for x in [1,2,3,4,5,6,7,8,9,10]]
    '[37.49, 60.22, 84.78, 110.97, 138.75, 168.16, 199.4, 232.92, 269.9, 296.89]'
    >>> [mity_qual(x,100) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[71.87, 156.08, 251.23, 354.34, 463.9, 579.05, 699.21, 824.04, 953.32, 1086.94]'
    >>> [mity_qual(x,1000) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[17.84, 50.97, 93.59, 142.89, 197.36, 256.02, 318.25, 383.57, 451.61, 522.1]'
    >>> [mity_qual(x,10000) for x in [5,10,15,20,25,30,35,40,45,50]]
    '[0.0, 0.05, 0.74, 3.56, 9.51, 18.73, 31.0, 46.04, 63.57, 83.37]'
    >>> mity_qual(118,118)
    '3211.76'
    >>> mity_qual(119,119)
    '3220'
    >>> mity_qual(10000,10000)
    '3220'
    """
    q = 0.0
    AO = int(AO)
    DP = int(DP)
    if AO > 0 and DP > 0:
        if DP == AO:
            # then penalise homoplasmic variants with low numbers of reads
            DP = DP + 1
        # v1: the cdf implementation caps at 159.55, due to the cdf capping at 0.9999999999999999
        # q = round(abs(-10 * log10(1 - binom.cdf(AO, DP, p))), 2)
        # v2: the logsf implementation is identical to cdf for low values and continues to scale up
        with np.errstate(divide='ignore'):
            q = float(round(-4.342945 * binom.logsf(AO, DP, p), 2))
    if isinf(q):
        # this is the maximum q that we observed looking at homoplastic variants around ~129x depth
        q = float(3220)
    return q

def add_filter(variant, allsamples, p):
    # filtering
    pass_flag = True
    if variant.pos in BLACKLIST:
        variant.filter.add("POS")
        pass_flag = False
    
    # filter dictionary
    # each value represents the number of samples that pass each field
    # starts with all samples passing
    num_samples = len(variant.samples.values())
    filter_dict = {
        "SBR": num_samples,
        "MQMR": num_samples,
        "AQR": num_samples,
        "SBA": num_samples
    }

    # testing samples
    for sample in variant.samples.values():
        # adding to format
        # AQR
        if int(sample["RO"]) > 0:
            AQR = float(sample["QR"]) / float(sample["RO"])
        else:
            AQR = 0
        
        sample["AQR"] = AQR

        # VAF
        if sample["DP"] != 0:
            VAF = round(float(sample["AO"][0]) / float(sample["DP"]), 4)
        else:
            VAF = 0
        
        sample["VAF"] = VAF

        # mity quality
        q = mity_qual(sample["AO"][0], sample["DP"], p)
        sample["q"] = q

        # AQA
        if sample["AO"][0] > 0:
            AQA = float(round(sample["QA"][0] / sample["AO"][0], 3))
        else:
            AQA = float(0)

        sample["AQA"] = AQA

        # tier (originally from mity report)
        if float(VAF) >= 0.01:
            tier = 1
        elif float(VAF) < 0.01 and float(sample["AO"][0]) > 10:
            tier = 2
        else:
            tier = 3

        sample["tier"] = tier

        # testing
        if sample["RO"] > MIN_DP:
            if not SB_RANGE_LO <= variant.info["SBR"] <= SB_RANGE_HI:
                filter_dict["SBR"] -= 1
            if variant.info["MQMR"] < MIN_MQMR:
                filter_dict["MQMR"] -= 1
            if AQR < MIN_AQR:
                filter_dict["AQR"] -= 1
        if sample["AO"][0] > MIN_DP:
            if not SB_RANGE_LO <= variant.info["SBA"] <= SB_RANGE_HI:
                filter_dict["SBA"] -= 1

    # all samples must pass each test
    if (allsamples):
        for (field, num_passing_samples) in filter_dict.items():
            if num_passing_samples != num_samples:
                variant.filter.add(field)
                pass_flag = False

    # only one sample has to pass each test
    else:
        for (field, num_passing_samples) in filter_dict.items():
            if num_passing_samples == 0:
                variant.filter.add(field)
                pass_flag = False

    # setting passing filter
    if pass_flag:
        variant.filter.add("PASS")

    return variant

def add_info_values(variant):
    # adding to INFO field

    # we want to caclulate SBR and add to INFO
    # SBR = strand bias of reference
    # SBR = SRF/(SRF+SRR)
    SRF = variant.info["SRF"]
    SRR = variant.info["SRR"]
    # to avoid dividing by zero:
    if SRF > 0 or SRR > 0:
        SBR = round(SRF / (SRF + SRR), 3)
    else:
        SBR = 0

    # we want to calculate SBA and add to INFO
    # SBA = strand bias of alternate
    # i.e. SBA = SAF/(SAF+SAR)

    # >>> variant.info["SAF"]
    # >>> (8505,)
    # so we use [0] to get the first value, and we can assume that there is only one value
    SAF = variant.info["SAF"][0]
    SAR = variant.info["SAR"][0]

    SBA = round(SAF / (SAF + SAR), 3)
        
    variant.info.update({
        "SBR": float(SBR), 
        "SBA": float(SBA)
    })

    return variant


def do_gsort(filtered_vcf, genome_file, output_file):
    gsort_cmd = "gsort {} {} | bgzip -cf > {}".format(filtered_vcf, genome_file, output_file)
    subprocess.run(gsort_cmd, shell=True)


def do_normalise(debug, vcf, reference_fasta, genome, output_file=None, allsamples=False, p=P_VAL):
    if debug:
        logger.setLevel(logging.DEBUG)
        logger.debug("Entered debug mode.")
    else:
        logger.setLevel(logging.INFO)

    normalised_vcf = vcf.replace(".vcf.gz", ".bcftoolsnorm.vcf.gz")

    # The "-o" option doesn't get passed to pysam's bcftools wrapper (as of 31-08-2023)
    # so we write to a separate file manually.
    with open(normalised_vcf, "w") as f:
        print(pysam.bcftools.norm("-f", reference_fasta, "-m-both", vcf), end="", file=f)
    
    # vcf files
    input_vcf = pysam.VariantFile(normalised_vcf)
    new_header = add_headers(input_vcf)
    filtered_vcf = vcf.replace(".vcf.gz", ".filtered.vcf")
    output_vcf = pysam.VariantFile(filtered_vcf, 'w', header=new_header)
    
    for variant in input_vcf.fetch():
        variant = add_info_values(variant)
        variant = add_filter(variant, allsamples, p)

        output_vcf.write(variant)
    output_vcf.close()
    input_vcf.close()

    # gsort
    if output_file is None:
        output_file = vcf.replace(".vcf.gz", ".mity.normalise.vcf.gz")

    do_gsort(filtered_vcf, genome, output_file)
    tabix(output_file)

    # remove temporary files
    os.remove(filtered_vcf)
    os.remove(normalised_vcf)