import pysam.bcftools

input_vcf_file = 'bcftools_normalise_test.vcf'
filtered_vcf_file = 'filtered.vcf'

filter_expression = '"POS=310"'

input_vcf = pysam.VariantFile(input_vcf_file)


# doesn't correctly exclude
filtered_file = pysam.bcftools.filter('--exclude', filter_expression, '--soft-filter', '\"POS\"', input_vcf_file)
print(filtered_file)

input_vcf.close()