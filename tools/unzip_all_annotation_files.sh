#!/bin/bash

# Make sure to run from inside the tools/ directory

DIR1="../mitylib/annot_chrm"
DIR2="../mitylib/annot_mt"

process_dir() {
    local DIR="$1"
    echo "Processing directory: $DIR"
    
    find "$DIR" -type f -name "*.vcf.gz" -print0 | while IFS= read -r -d '' file; do
        echo "Decompressing: $file"
        gunzip -k "$file"
    done
}

process_dir "$DIR1"
process_dir "$DIR2"