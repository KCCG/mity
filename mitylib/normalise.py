"""Normalise mity VCF."""

import logging
import os.path
from math import isinf
from scipy.stats import binom
import numpy as np
import pysam
import pysam.bcftools
from mitylib.util import MityUtil


# LOGGING
logger = logging.getLogger(__name__)


class Normalise:
    """
    Mity normalise.
    """

    P_VAL = 0.002
    SB_RANGE_LO = 0.1
    SB_RANGE_HI = 0.9
    MIN_MQMR = 30
    MIN_AQR = 20
    MIN_DP = 15
    BLACKLIST = list(range(302, 319)) + list(range(3105, 3109))

    def __init__(
        self,
        debug,
        vcf,
        reference_fasta,
        genome,
        output_dir=".",
        prefix=None,
        allsamples=False,
        keep=False,
        p=P_VAL,
    ):
        self.debug = debug
        self.vcf = vcf
        self.reference_fasta = reference_fasta
        self.genome = genome
        self.output_dir = output_dir
        self.prefix = prefix
        self.allsamples = allsamples
        self.keep = keep
        self.p = p

        self.bcftools_norm_obj = None
        self.filtered_vcf_obj = None

        self.bcftools_norm_path = ""
        self.filtered_vcf_path = ""
        self.normalised_vcf_path = ""

        self.run()

    def run(self):
        """
        Run mity normalise.
        """
        if self.debug:
            logger.setLevel(logging.DEBUG)
            logger.debug("Entered debug mode.")
        else:
            logger.setLevel(logging.INFO)

        self.set_strings()
        self.run_bcftools_norm()
        self.run_filtering()

        MityUtil.gsort(self.filtered_vcf_path, self.normalised_vcf_path, self.genome)

        self.remove_intermediate_files()

    def run_filtering(self):
        """
        Makes pysam VariantFile objects and runs filtering for each variant.
        """
        self.bcftools_norm_obj = pysam.VariantFile(self.bcftools_norm_path)
        self.filtered_vcf_obj = pysam.VariantFile(
            self.filtered_vcf_path, "w", header=self.add_headers()
        )

        for variant in self.bcftools_norm_obj.fetch():
            variant = self.add_info_values(variant)
            variant = self.add_filter(variant)

            self.filtered_vcf_obj.write(variant)

        self.bcftools_norm_obj.close()
        self.filtered_vcf_obj.close()

    def run_bcftools_norm(self):
        """
        Run bcftools norm.
        """
        # The "-o" option doesn't get passed to pysam's bcftools wrapper (as of 31-08-2023)
        # so we write to a separate file manually.
        with open(self.bcftools_norm_path, "w") as f:
            print(
                pysam.bcftools.norm("-f", self.reference_fasta, "-m-both", self.vcf),
                end="",
                file=f,
            )

    def set_strings(self):
        """
        Sets:
            prefix
            normalised_vcf
            filtered_vcf
        """
        if self.prefix is None:
            self.prefix = MityUtil.make_prefix(self.vcf)

        self.bcftools_norm_path = os.path.join(
            self.output_dir, self.prefix + ".bcftools.norm.vcf.gz"
        )
        self.filtered_vcf_path = self.vcf.replace(".vcf.gz", ".filtered.vcf")
        self.normalised_vcf_path = os.path.join(
            self.output_dir, self.prefix + ".normalise.vcf.gz"
        )

    def add_headers(self):
        """
        Return new headers for mity normalise vcf output.
        """
        new_header = self.bcftools_norm_obj.header

        # fitler headers
        new_header.filters.add(
            "FAIL",
            None,
            None,
            "Variant fails filter in at least one or all samples (depending on filter option)",
        )

        # info headers
        new_header.info.add(
            "SBR",
            number=1,
            type="Float",
            description="For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)",
        )
        new_header.info.add(
            "SBA",
            number=1,
            type="Float",
            description="For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)",
        )

        # format headers
        new_header.formats.add(
            "AQA",
            number="A",
            type="Float",
            description="Average base quality of the alternate reads, AQA=QA/AO",
        )
        new_header.formats.add(
            "AQR",
            number="A",
            type="Float",
            description="Average base quality of the reference reads, AQR=QR/RO",
        )
        new_header.formats.add(
            "VAF",
            number="A",
            type="Float",
            description="Allele frequency in the range (0,1] - the ratio of the number of alternate reads to reference reads",
        )
        new_header.formats.add(
            "q",
            number="A",
            type="Float",
            description="Phred scaled binomial probability of seeing AO reads from DP, assuming a noise floor of p=0.002. ",
        )
        new_header.formats.add(
            "tier",
            number="A",
            type="Integer",
            description="Custom variant tier (TODO: change description)",
        )

        # format filter headers
        new_header.formats.add(
            "POS_filter",
            number="A",
            type="Integer",
            description="Variant falls in the blacklist of positions: MT:302-319, MT:3105-3108",
        )
        new_header.formats.add(
            "SBR_filter",
            number="A",
            type="Integer",
            description="For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)",
        )
        new_header.formats.add(
            "SBA_filter",
            number="A",
            type="Integer",
            description="For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)",
        )
        new_header.formats.add(
            "MQMR_filter",
            number="A",
            type="Integer",
            description="For all alleles MQMR<30",
        )
        new_header.formats.add(
            "AQR_filter",
            number="A",
            type="Integer",
            description="For all alleles AQR<20",
        )

        return new_header

    def mity_qual(self, AO, DP, p):
        """
        Calculate the Phred-scaled quality score for a genetic variant.

        Parameters:
            - AO (int): Number of alternative reads.
            - DP (int): Total read depth.

        Returns:
            - float: Phred-scaled quality score.
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
            with np.errstate(divide="ignore"):
                q = float(round(-4.342945 * binom.logsf(AO, DP, p), 2))
        if isinf(q):
            # this is the maximum q that we observed looking at homoplastic variants around ~129x depth
            q = float(3220)
        return q

    def add_filter(self, variant):
        """
        Adds filter to variant, and filters for samples.

        Parameters:
            - variant (VariantRecord): Variant to add filters to

        Returns:
            - VariantRecord: Variant with filters added
        """
        # filter dictionary
        # each value represents the number of samples that pass each field
        # starts with all samples passing
        num_samples = len(variant.samples.values())
        filter_dict = {
            "SBR": num_samples,
            "MQMR": num_samples,
            "AQR": num_samples,
            "SBA": num_samples,
        }

        # filtering
        pass_flag = True
        pos_flag = True
        if variant.pos in self.BLACKLIST:
            pos_flag = False
            pass_flag = False

        # filtering samples
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
            q = self.mity_qual(sample["AO"][0], sample["DP"], self.p)
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

            # set tests to "PASS", i.e. 1 first
            # sidenote: it would be nice to use TRUE/FALSE but the format field does
            # not support boolean values
            sample["POS_filter"] = 1 if pos_flag else 0
            sample["SBR_filter"] = 1
            sample["SBA_filter"] = 1
            sample["MQMR_filter"] = 1
            sample["AQR_filter"] = 1

            # testing
            if sample["RO"] > self.MIN_DP:
                if not self.SB_RANGE_LO <= variant.info["SBR"] <= self.SB_RANGE_HI:
                    filter_dict["SBR"] -= 1
                    sample["SBR_filter"] = 0
                if variant.info["MQMR"] < self.MIN_MQMR:
                    filter_dict["MQMR"] -= 1
                    sample["MQMR_filter"] = 0
                if AQR < self.MIN_AQR:
                    filter_dict["AQR"] -= 1
                    sample["AQR_filter"] = 0
            if sample["AO"][0] > self.MIN_DP:
                if not self.SB_RANGE_LO <= variant.info["SBA"] <= self.SB_RANGE_HI:
                    filter_dict["SBA"] -= 1
                    sample["SBA_filter"] = 0

        # all samples must pass each test
        if self.allsamples:
            for num_passing_samples in filter_dict.values():
                if num_passing_samples != num_samples:
                    pass_flag = False
                    break

        # only one sample has to pass each test
        else:
            for num_passing_samples in filter_dict.values():
                if num_passing_samples == 0:
                    pass_flag = False
                    break

        # setting passing filter
        if pass_flag:
            variant.filter.add("PASS")
        else:
            variant.filter.add("FAIL")

        return variant

    def add_info_values(self, variant):
        """
        Adds new and calculated INFO values to a variant.
        """
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

        variant.info.update({"SBR": float(SBR), "SBA": float(SBA)})

        return variant

    def remove_intermediate_files(self):
        """
        Remove intermediate files.
        """
        if not self.keep:
            os.remove(self.filtered_vcf_path)
            os.remove(self.bcftools_norm_path)
