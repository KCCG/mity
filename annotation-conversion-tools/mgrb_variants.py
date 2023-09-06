def convert_to_vcf(input_file, output_file):
    vcf_header = [
        '##fileformat=VCFv4.2',
        '##INFO=<ID=MGRB_FILTER,Number=1,Type=String,Description="MGRB Filter">',
        '##INFO=<ID=MGRB_AN,Number=1,Type=Integer,Description="MGRB AN">',
        '##INFO=<ID=MGRB_AC,Number=A,Type=Integer,Description="MGRB AC">',
        '##INFO=<ID=MGRB_frequency,Number=A,Type=Float,Description="MGRB Frequency">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]

    vcf_data = []
    with open(input_file, 'r') as infile:
        for line in infile:
            fields = line.strip().split(',')
            chrom = fields[0]
            pos = fields[1]
            ref = fields[2]
            alt = fields[3]
            mgrb_filter = fields[4]
            mgrb_an = fields[5]
            mgrb_ac = fields[6]
            mgrb_frequency = fields[7]

            vcf_data.append(f'{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\tMGRB_FILTER={mgrb_filter};MGRB_AN={mgrb_an};MGRB_AC={mgrb_ac};MGRB_frequency={mgrb_frequency}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(vcf_header + vcf_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/mgrb_variants.csv"  # Change this to the path of your input text file
    output_file = "output.vcf"  # Change this to the desired output VCF file path

    convert_to_vcf(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
