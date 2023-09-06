def convert_to_bed(input_file, output_file):
    bed_data = []

    with open(input_file, 'r') as infile:
        next(infile)  # Skip the header line
        for line in infile:
            fields = line.strip().split(',')
            chrom = fields[0]
            map_locus = fields[1]
            start_pos = int(fields[2]) - 1  # Convert 1-based position to 0-based for BED format
            end_pos = int(fields[2])  # End position is the same as the start position
            description = fields[3]

            bed_data.append(f'{chrom}\t{start_pos}\t{end_pos}\t{map_locus}\t.\t.\t{description}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(bed_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/mito_dna_func_loc.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/mito_dna_func_loc.bed"  # Change this to the desired output BED file path

    convert_to_bed(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
