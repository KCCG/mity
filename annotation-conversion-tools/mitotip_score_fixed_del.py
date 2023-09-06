import csv

def convert_to_vcf(input_file, output_file):
    # VCF header lines
    vcf_header = [
        '##fileformat=VCFv4.3',
        '##INFO=<ID=MitoTip_score,Number=1,Type=Float,Description="MitoTip Score">',
        '##INFO=<ID=MitoTip_percentile,Number=1,Type=Float,Description="MitoTip Percentile">',
        '##INFO=<ID=MitoTip_interpretation,Number=1,Type=String,Description="MitoTip Interpretation">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]

    vcf_data = []
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        for row in reader:
            pos = int(row['POS'])
            ref = row['REF']
            alt = row['ALT']
            score = float(row['MitoTip_score'])
            percentile = float(row['MitoTip_percentile'])
            interpretation = row['MitoTip_interpretation']

            info_fields = [
                f'MitoTip_score={score}',
                f'MitoTip_percentile={percentile}',
                f'MitoTip_interpretation={interpretation}'
            ]

            info = ';'.join(info_fields)

            vcf_data.append(f'MT\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(vcf_header + vcf_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/mitotip_score_fixed_del.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/mitotip_score_fixed_del.vcf"  # Change this to the desired output VCF file path

    convert_to_vcf(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
