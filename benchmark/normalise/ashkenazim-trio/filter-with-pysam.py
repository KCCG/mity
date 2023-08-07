import pysam

# CONSTANTS
P_VAL = 0.002
SB_RANGE_LO = 0.1 
SB_RANGE_HI = 0.9
MIN_MQMR = 30
MIN_AQR = 20
MIN_DP = 15

BLACKLIST = list(range(302, 319)) + list(range(3105,3109))

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
    new_header.info.add("SBR", number=1, type="Integer", description="For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)")
    new_header.info.add("SBA", number=1, type="Integer", description="For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)")
    new_header.info.add("AQR", number=1, type="Integer", description="For all alleles AQR<20")

    return new_header


def add_filter(input_vcf, min_DP=MIN_DP, SB_range_lo=SB_RANGE_LO, SB_range_hi=SB_RANGE_HI, min_MQMR=MIN_MQMR, min_AQR=MIN_AQR):
    # setting headers
    new_header = input_vcf.header

    # filter headers
    new_header.filters.add("POS", None, None, "Variant falls in the blacklist of positions: MT:302-319, MT:3105-3108")
    new_header.filters.add("SBR", None, None, "For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)")
    new_header.filters.add("SBA", None, None, "For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)")
    new_header.filters.add("MQMR", None, None, "For all alleles MQMR<30")
    new_header.filters.add("AQR", None, None, "For all alleles AQR<20")

    # info headers
    new_header.info.add("SBR", number=1, type="Integer", description="For all alleles RO > 15 and (SBR > 0.9 or SBR < 0.1)")
    new_header.info.add("SBA", number=1, type="Integer", description="For all alleles AO > 15 and (SBA > 0.9 or SBA < 0.1)")
    new_header.info.add("AQR", number=1, type="Integer", description="For all alleles AQR<20")

    # output vcf file
    output_vcf = pysam.VariantFile('filtered.vcf', 'w', header=new_header)
    
    for variant in input_vcf.fetch():
        # TODO: move this somewhere else, potentially use vcfanno to annotate before filtering
        # adding SBR and SBA to INFO field

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

        # TODO: figure out a way to add something to the info column
        # variant.info.update({"SBR": SBR})

        # we want to calculate SBA and add to INFO
        # SBA = strand bias of alternate
        # i.e. SBA = SAF/(SAF+SAR)

        # >>> variant.info["SAF"]
        # >>> (8505,)
        # so we use [0] to get the first value, and we can assume that there is only one value
        SAF = variant.info["SAF"][0]
        SAR = variant.info["SAR"][0]

        SBA = round(SAF / (SAF + SAR), 3)
        # variant.info.update({"SBA": SBA})

        # filtering
        if variant.pos in BLACKLIST:
            variant.filter.add("POS")
        
        # currently requires all samples to pass the filter, (different from mity normalise)
        # TODO: add options for only requiring one sample to pass vs all vs percentage
        for sample in variant.samples.values():
            # AQR: note this is supposed to be added to the format + samples section
            if int(sample["RO"]) > 0:
                AQR = float(sample["QR"]) / float(sample["RO"])
            else:
                AQR = 0

            if sample["RO"] > min_DP:
                if not SB_range_lo <= SBR <= SB_range_hi:
                    variant.filter.add("SBR")
                if variant.info["MQMR"] < min_MQMR:
                    variant.filter.add("MQMR")
                if AQR < min_AQR:
                    variant.filter.add("AQR")
            if sample["AO"][0] > min_DP:
                if not SB_range_lo <= SBA <= SB_range_hi:
                    variant.filter.add("SBA")

        # if no other filters have been added to filter, then the variant passes
        if len(variant.filter.values()) == 0:
            variant.filter.add("PASS")

        output_vcf.write(variant)
        
    output_vcf.close()


def do_normalise():
    input_vcf = pysam.VariantFile('bcftools_normalise_test.vcf')
    # TODO: make a class/struct for the values for minDP, SBR, etc? 
    add_filter(input_vcf, min_DP=MIN_DP, SB_range_lo=SB_RANGE_LO, SB_range_hi=SB_RANGE_HI, min_MQMR=MIN_MQMR, min_AQR=MIN_AQR)
    input_vcf.close()

def main():
    do_normalise()

if __name__ == '__main__':
    main()