#!/bin/bash

# Check if the folder path is provided
if [ -z "$1" ]; then
    echo "Usage: $0 <folder_path>"
    exit 1
fi

FOLDER_PATH=$1

# Remove existing .gz and .tbi files
find "$FOLDER_PATH" -type f \( -name "*.bed.gz.tbi" -o -name "*.bed.gz" -o -name "*.vcf.gz.tbi" -o -name "*.vcf.gz" \) -exec rm -f {} +

# Bgzip and tabix .vcf and .bed files
for file in "$FOLDER_PATH"/*.{vcf,bed}; do
    if [ -f "$file" ]; then
        bgzip -k "$file"
        tabix "$file.gz"
    fi
done