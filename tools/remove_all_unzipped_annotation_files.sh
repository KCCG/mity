#!/bin/bash

# Make sure to run from inside the tools/ directory

DIR1="../mitylib/annot_chrm"
DIR2="../mitylib/annot_mt"

process_dir() {
    local DIR="$1"
    echo "Processing directory: $DIR"

    find "$DIR" -type f \( -name "*.vcf" -o -name "*.bed" \) -print0 | while IFS= read -r -d '' file; do
        echo "Removing: $file"
        rm "$file"
    done
}

process_dir "$DIR1"
process_dir "$DIR2"
