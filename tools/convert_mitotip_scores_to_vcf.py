"""
Script to convert mitotip scores to VCF format.

From: https://www.mitomap.org/MITOMAP/MitoTipScores

Usage: python convert_mitotip_scores_to_vcf.py mitomap_scores.txt mitomap_scores.vcf
"""

import sys
import csv
import pysam


def quartile_to_percentile(quartile):
    return {"Q1": "0.75-1.00", "Q2": "0.50-0.75", "Q3": "0.25-0.50", "Q4": "0.00-0.25"}.get(
        quartile, "NOT_PRESENT"
    )


def quartile_to_interpretation(quartile):
    return {
        "Q1": "likely-pathogenic",
        "Q2": "possibly-pathogenic",
        "Q3": "possibly-benign",
        "Q4": "likely-benign",
    }.get(quartile, "NOT_PRESENT")


def parse_tsv_to_vcf(input_tsv, output_vcf):
    # Define VCF header
    vcf_header = pysam.VariantHeader()
    vcf_header.add_line("##fileformat=VCFv4.3")

    # Contig
    # Change between MT and chrM as needed
    # Remember to change this below too!
    vcf_header.add_line("##contig=<ID=MT,length=16569>")
    # vcf_header.add_line("##contig=<ID=chrM,length=16569>")

    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "MitoTip_score"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description", "MitoTip Score"),
        ],
    )
    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "MitoTip_percentile"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "MitoTip Percentile"),
        ],
    )
    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "MitoTip_quartile"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "MitoTip Quartile"),
        ],
    )
    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "MitoTip_score_interpretation"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "MitoTip Score Interpretation"),
        ],
    )

    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "Mitomap_status"),
            ("Number", "1"),
            ("Type", "String"),
            ("Description", "Mitomap Status"),
        ],
    )
    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "Count"),
            ("Number", "1"),
            ("Type", "Integer"),
            ("Description", "Observation Count"),
        ],
    )
    vcf_header.add_meta(
        "INFO",
        items=[
            ("ID", "Percentage"),
            ("Number", "1"),
            ("Type", "Float"),
            ("Description", "Observation Percentage"),
        ],
    )
    vcf_header.add_sample("SAMPLE")  # VCF requires at least one sample, even if not used

    # Open VCF file for writing
    with (
        pysam.VariantFile(output_vcf, "w", header=vcf_header) as vcf_out,
        open(input_tsv, "r", encoding="utf-8") as infile,
    ):
        reader = csv.DictReader(infile, delimiter="\t")

        for row in reader:
            # Change between MT and chrM as needed
            chrom = "MT"
            # chrom = "chrM"

            pos = int(row["Position"])
            ref = row["rCRS"]
            alt = row["Alt"]
            score = float(row["MitoTIP_Score"])
            quartile = row["Quartile"]
            count = int(row["Count"])
            percentage = float(row["Percentage"])
            status = row["Mitomap_Status"]

            # Convert quartile to percentile
            percentile = quartile_to_percentile(quartile)

            # Convert quartile to interpretation
            interpretation = quartile_to_interpretation(quartile)

            # Create new record
            record = vcf_out.new_record()

            # Set fields manually
            record.chrom = chrom
            record.pos = pos
            record.id = "."  # No specific ID
            record.ref = ref
            record.alts = (alt,)  # Must be a tuple
            record.qual = None  # No quality score
            record.filter.add("PASS")  # Default to PASS filter
            record.info["MitoTip_score"] = score
            record.info["MitoTip_percentile"] = percentile
            record.info["MitoTip_quartile"] = quartile
            record.info["MitoTip_score_interpretation"] = interpretation
            record.info["Mitomap_status"] = status
            record.info["Count"] = count
            record.info["Percentage"] = percentage

            # Remove FILTER field (make it empty)
            record.filter.clear()

            # Write to VCF
            vcf_out.write(record)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python convert_to_vcf.py input.tsv output.vcf")
        sys.exit(1)

    input_tsv = sys.argv[1]
    output_vcf = sys.argv[2]

    parse_tsv_to_vcf(input_tsv, output_vcf)
