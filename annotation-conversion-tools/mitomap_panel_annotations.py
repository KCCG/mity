def convert_to_vcf(input_file, output_file):
    vcf_header = [
        '##fileformat=VCFv4.3',
        '##INFO=<ID=locus_mitomap,Number=1,Type=String,Description="Locus from Mitomap">',
        '##INFO=<ID=num_references_mitomap,Number=1,Type=Integer,Description="Number of references from Mitomap">',
        '##INFO=<ID=variant_amino_acid_change_mitomap,Number=1,Type=String,Description="Variant amino acid change from Mitomap">',
        '##INFO=<ID=codon_position_mitomap,Number=1,Type=String,Description="Codon position from Mitomap">',
        '##INFO=<ID=codon_number_mitomap,Number=1,Type=String,Description="Codon number from Mitomap">',
        '##INFO=<ID=disease_mitomap,Number=1,Type=String,Description="Disease from Mitomap">',
        '##INFO=<ID=num_disease_references_mitomap,Number=1,Type=Integer,Description="Number of disease references from Mitomap">',
        '##INFO=<ID=RNA_mitomap,Number=1,Type=String,Description="RNA from Mitomap">',
        '##INFO=<ID=homoplasmy_mitomap,Number=1,Type=String,Description="Homoplasmy from Mitomap">',
        '##INFO=<ID=heteroplasmy_mitomap,Number=1,Type=String,Description="Heteroplasmy from Mitomap">',
        '##INFO=<ID=status_mitomap,Number=1,Type=String,Description="Status from Mitomap">',
        '##INFO=<ID=disease_amino_acid_change_mitomap,Number=1,Type=String,Description="Disease amino acid change from Mitomap">',
        '##INFO=<ID=GenBank_frequency_mitomap,Number=1,Type=String,Description="GenBank frequency from Mitomap">',
        '##INFO=<ID=baylor_panel,Number=1,Type=String,Description="Baylor panel from Mitomap">',
        '##INFO=<ID=common_22_panel,Number=1,Type=String,Description="Common 22 panel from Mitomap">',
        '##INFO=<ID=common_58_panel,Number=1,Type=String,Description="Common 58 panel from Mitomap">',
        '##INFO=<ID=commercial_panels,Number=1,Type=String,Description="Commercial panels from Mitomap">',
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO'
    ]

    vcf_data = []

    with open(input_file, 'r') as infile:
        next(infile)  # Skip the header line
        for line in infile:
            fields = line.strip().split(',')
            pos = int(fields[0])
            ref = fields[1]
            alt = fields[2]
            info_fields = [
                f'locus_mitomap={fields[3]}',
                f'num_references_mitomap={fields[4]}',
                f'variant_amino_acid_change_mitomap={fields[5]}',
                f'codon_position_mitomap={fields[6]}',
                f'codon_number_mitomap={fields[7]}',
                f'disease_mitomap={fields[8]}',
                f'num_disease_references_mitomap={fields[9]}',
                f'RNA_mitomap={fields[10]}',
                f'homoplasmy_mitomap={fields[11]}',
                f'heteroplasmy_mitomap={fields[12]}',
                f'status_mitomap={fields[13]}',
                f'disease_amino_acid_change_mitomap={fields[14]}',
                f'GenBank_frequency_mitomap={fields[15]}',
                f'baylor_panel={fields[16]}',
                f'common_22_panel={fields[17]}',
                f'common_58_panel={fields[18]}',
                f'commercial_panels={fields[20]}'
            ]

            info = ';'.join(info_fields)

            vcf_data.append(f'MT\t{pos}\t.\t{ref}\t{alt}\t.\t.\t{info}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(vcf_header + vcf_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/mitomap_panel_annotations_sorted.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/mitomap_panel_annotations.vcf"  # Change this to the desired output VCF file path

    convert_to_vcf(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
