import csv
import sys
import pysam
import re


def parse_allele(allele):
    match = re.match(r"m\.(\d+)([A-Zdel]+)>([A-Z]*)", allele)
    if match:
        pos, ref, alt = match.groups()
        return int(pos), ref, alt or "."
    return None, ".", "."

def remove_semicolon(value):
    return value.replace(";", ",")


def tsv_to_vcf(tsv_file, vcf_file):
    header = pysam.VariantHeader()
    header.add_line("##fileformat=VCFv4.2")
    header.add_line("##contig=<ID=MT,length=16569>")

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
            ("ID", "LocusType"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Locus Type"),
        ],
    )

    header.add_meta(
        "INFO",
        items=[
            ("ID", "AssociatedDisease"),
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
        items=[
            ("ID", "AminoAcidChange"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Amino Acid Change"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "StatusMitomapClinGen"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Status Mitomap ClinGen"),
        ],
    )
    header.add_meta(
        "INFO",
        items=[
            ("ID", "LastUpdate"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Last Update"),
        ],
    )

    header.add_sample("SAMPLE")

    with (
        pysam.VariantFile(vcf_file, "w", header=header) as vcfout,
        open(tsv_file, newline="", encoding="utf-8") as tsvfile,
    ):
        reader = csv.DictReader(tsvfile, delimiter=",")

        for row in reader:
            pos, ref, alt = parse_allele(row["Allele"])
            if pos is None:
                continue

            record = vcfout.new_record()
            record.contig = "MT"
            record.start = pos - 1
            record.stop = pos
            record.id = "."
            record.ref = ref
            record.alts = [alt] if alt != "." else None
            record.qual = None
            record.filter.add(".")

            info_fields = {
                "Locus": row["Locus"],
                "LocusType": row["Locus Type"],
                "AssociatedDisease": row["Associated Diseases"],
                "Allele": row["Allele"],
                "AminoAcidChange": row["aaΔ\xa0or\xa0RNA"],
                "StatusMitomapClinGen": row["Status\xa0♣(Mitomap [ClinGen])"],
                "LastUpdate": row["Last\xa0StatusUpdate"]
            }

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
