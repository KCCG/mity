import pysam

input_vcf = pysam.VariantFile('bcftools_normalise_test.vcf')

# setting headers
new_header = input_vcf.header
new_header.filters.add("FAIL", None, None, "description")
new_header.filters.add("TEST", None, None, "description")

output_vcf = pysam.VariantFile('filtered.vcf', 'w', header=new_header)

BLACKLIST = list(range(302, 319)) + list(range(3105,3109))

for variant in input_vcf.fetch():
    # POS
    # if variant.pos in BLACKLIST:
    #     variant.filter.add("FAIL")
    # else:
    #     variant.filter.add("PASS")
    variant.filter.add("FAIL")
    variant.filter.add("TEST")

    output_vcf.write(variant)

input_vcf.close()
output_vcf.close()