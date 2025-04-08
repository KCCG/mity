import csv
import sys
import pysam


def parse_allele(allele):
    """Extract REF and ALT from allele notation like A576G."""
    if len(allele) < 2:
        return ".", "."
    ref = allele[0]  # First character is REF
    alt = allele[-1]  # Last character is ALT
    return ref, alt


def remove_semicolon(value):
    return value.replace(";", "")


def csv_to_vcf(csv_file, vcf_file):
    """Converts a CSV file into a VCF file using pysam."""
    header = pysam.VariantHeader()

    # Standard VCF fields
    header.add_line("##fileformat=VCFv4.2")
    # Contig
    # Change between MT and chrM as needed
    # Remember to change this below too!
    header.add_line("##contig=<ID=MT,length=16569>")
    # vcf_header.add_line("##contig=<ID=chrM,length=16569>")

    # Define INFO fields for all columns
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Locus"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Locus"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Disease"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Associated Disease"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Allele"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Allele"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[("ID", "RNA"), ("Number", "1"), ("Type", "String"), ("Description", "RNA Type")],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Homoplasmy"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Homoplasmy Status"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Heteroplasmy"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Heteroplasmy Status"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "Status"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Clinical Status"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "MitoTIP"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "MitoTIP Score"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "GB_Freq"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "GenBank Frequency"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "GB_Seqs"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "GenBank Sequences"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "References"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Reference Count"),
        ],
    )

    header.add_sample("SAMPLE")

    # Create VCF file
    with (
        pysam.VariantFile(vcf_file, "w", header=header) as vcfout,
        open(csv_file, newline="", encoding="utf-8") as csvfile,
    ):
        reader = csv.DictReader(csvfile)

        for row in reader:
            chrom = "MT"  # Mitochondrial chromosome
            pos = int(row["Position"])
            id_ = "."
            ref, alt = parse_allele(row["Allele"])
            qual = None  # No quality scores

            # Construct INFO fields
            info_fields = {
                "Locus": row["Locus"],
                "Disease": row["Disease"],
                "Allele": row["Allele"],
                "RNA": row["RNA"],
                "Homoplasmy": row["Homoplasmy"],
                "Heteroplasmy": row["Heteroplasmy"],
                "Status": row["Status"],
                "MitoTIP": row["MitoTIP†"],
                "GB_Freq": row["GB\xa0Freq\xa0\xa0FL\xa0(CR)*‡"],
                "GB_Seqs": row["GB\xa0Seqs\xa0FL\xa0(CR)*"],
                "References": row["References"],
            }

            # Create a new VCF record
            record = vcfout.new_record()
            record.contig = chrom
            record.start = pos - 1  # VCF uses 0-based positions
            record.stop = pos
            record.id = id_
            record.ref = ref
            record.alts = [alt] if alt != "." else None
            record.qual = qual
            record.filter.add(".")

            for key, value in info_fields.items():
                record.info[key] = remove_semicolon(value)

            # Write the record
            vcfout.write(record)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_vcf.py input.tsv output.vcf")
        sys.exit(1)

    input_csv = sys.argv[1]
    output_vcf = sys.argv[2]

    csv_to_vcf(input_csv, output_vcf)
