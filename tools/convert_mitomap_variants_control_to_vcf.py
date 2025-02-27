import csv
import sys
import pysam


def remove_semicolon(value):
    return value.replace(";", ",")


def parse_nucleotide_change(change):
    ref, alt = change.split("-")
    if alt.lower() == "del":
        return ref, "."  # Deletion
    return ref, alt  # Substitution or insertion


def tsv_to_vcf(tsv_file, vcf_file):
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.add_line("##contig=<ID=MT,length=16569>")

    # Define INFO fields based on the CSV columns
    info_field_headings = {
        "Locus": ("1", "String", "Locus"),
        "NucleotideChange": ("1", "String", "Nucleotide Change"),
        "CodonNumber": ("1", "String", "Codon Number"),
        "CodonPosition": ("1", "String", "Codon Position"),
        "AminoAcidChange": ("1", "String", "Amino Acid Change"),
        "GB_Freq": ("1", "String", "GenBank Frequency"),
        "GB_Seqs": ("1", "String", "GenBank Sequences"),
        "CuratedReferences": ("1", "String", "Curated References"),
    }

    for key, (num, typ, desc) in info_field_headings.items():
        header.add_meta(
            "INFO", items=[("ID", key), ("Number", num), ("Type", typ), ("Description", desc)]
        )

    header.add_sample("SAMPLE")

    with (
        pysam.VariantFile(vcf_file, "w", header=header) as vcfout,
        open(tsv_file, newline="", encoding="utf-8") as tsvfile,
    ):
        reader = csv.DictReader(tsvfile)

        for row in reader:
            pos = int(row["Position"])  # Convert Position to integer
            ref, alt = parse_nucleotide_change(row["Nucleotide Change"])
            
            print(ref, alt)

            record = vcfout.new_record()
            record.contig = "MT"
            record.start = pos - 1
            record.stop = pos
            record.id = "."
            record.ref = ref
            record.alts = [alt]
            record.qual = None
            record.filter.add(".")

            info_fields = {
                "Locus": row["Locus"],
                "NucleotideChange": row["Nucleotide Change"],
                "CodonNumber": row["Codon Number"],
                "CodonPosition": row["Codon Position"],
                "AminoAcidChange": row["Amino Acid Change"],
                "GB_Freq": row["GB Freq‡"],
                "GB_Seqs": row["GB Seqs"],
                "CuratedReferences": row["Curated References"],
            }

            # Fill INFO fields
            for key, value in info_fields.items():
                record.info[key] = remove_semicolon(value)

            vcfout.write(record)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_vcf.py input.tsv output.vcf")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_vcf = sys.argv[2]
    tsv_to_vcf(input_tsv, output_vcf)
