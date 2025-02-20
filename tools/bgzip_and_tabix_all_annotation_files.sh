#!/bin/bash

# Make sure to run from inside the tools/ directory

DIR1="../mitylib/annot_chrm"
DIR2="../mitylib/annot_mt"

process_dir() {
    local DIR="$1"
    echo "Processing directory: $DIR"

    find "$DIR" -type f -name "*.vcf.gz" -delete
    find "$DIR" -type f -name "*.vcf.gz.tbi" -delete
    find "$DIR" -type f -name "*.bed.gz" -delete
    find "$DIR" -type f -name "*.bed.gz.tbi" -delete

    find "$DIR" -type f -name "*.vcf" -print0 | while IFS= read -r -d '' file; do
        bgzip -k "$file"
        tabix -p vcf "${file}.gz"
    done

    find "$DIR" -type f -name "*.bed" -print0 | while IFS= read -r -d '' file; do
        bgzip -kf "$file"
        tabix -p bed -f "${file}.gz"
    done
}

process_dir "$DIR1"
process_dir "$DIR2"
