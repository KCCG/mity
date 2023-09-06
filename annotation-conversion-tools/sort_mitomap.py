import csv

def sort_csv_by_position(input_file, output_file):
    # Read the CSV file and sort it by the "POS" column
    with open(input_file, 'r') as infile:
        reader = csv.DictReader(infile)
        sorted_rows = sorted(reader, key=lambda row: int(row['POS']))

    # Write the sorted data to a new CSV file
    with open(output_file, 'w', newline='') as outfile:
        fieldnames = reader.fieldnames
        writer = csv.DictWriter(outfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(sorted_rows)

if __name__ == "__main__":
    input_file = "../mitylib/annot/mitomap_panel_annotations.csv"  # Change this to the path of your input CSV file
    output_file = "../mitylib/annot/mitomap_panel_annotations_sorted.csv"  # Change this to the desired output sorted CSV file path

    sort_csv_by_position(input_file, output_file)
    print(f"Sorting completed. Sorted output written to {output_file}")
