def convert_to_bed(input_file, output_file):
    bed_data = []

    with open(input_file, 'r') as infile:
        next(infile)  # Skip the header line
        for line in infile:
            fields = line.strip().split(',')
            chrom = fields[0]
            start_pos = int(fields[1]) - 1  # Convert 1-based position to 0-based for BED format
            end_pos = int(fields[1])  # End position is the same as the start position
            gene_name = fields[2]
            gene_biotype = fields[3]

            bed_data.append(f'{chrom}\t{start_pos}\t{end_pos}\t{gene_name}\t{gene_biotype}')

    with open(output_file, 'w') as outfile:
        outfile.write('\n'.join(bed_data))

if __name__ == "__main__":
    input_file = "../mitylib/annot/gtf_annotations.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/gtf_annotations.bed"  # Change this to the desired output BED file path

    convert_to_bed(input_file, output_file)
    print(f"Conversion completed. Output written to {output_file}")
