def convert_to_vcf(input_file, output_file):
    # VCF header lines
    vcf_header = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=phylotree_mut,Number=1,Type=String,Description="Phylotree Mutation">',
        '##INFO=<ID=phylotree_haplotype,Number=1,Type=String,Description="Phylotree Haplotype">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]

    vcf_data = []
    with open(input_file, 'r') as infile:
        for line in infile:
            fields = line.strip().split(',')
            chrom = fields[-1]  # Assuming CHR is the last column
            pos = fields[3]
            ref = fields[1]
            alt = fields[2]
            phylotree_mut = fields[0]
            phylotree_haplotype = fields[4]

            vcf_data.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tphylotree_mut={phylotree_mut};phylotree_haplotype={phylotree_haplotype}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(vcf_header + vcf_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/haplotype_data.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/haplotype_data.vcf"  # Change this to the desired output VCF file path

    convert_to_vcf(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
