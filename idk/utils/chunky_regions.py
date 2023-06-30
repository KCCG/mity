#!/usr/bin/env python3

# This code is adapted from fasta_generate_regions.py from https://github.com/ekg/freebayes

# usage: python3 chunky_regions.py --chunk_size INT --region [CHR:START-END] --bam_header_path [PATH]"
# or 
# Usage: bam_header stdin | python3 chunky_regions.py --chunk_size INT --region [CHR[:START-END]]

# if no chunk_size provided, then a file with a single column of the chromsomes is produced

# to get the bam header:
# samtools view -H <bam_path>

import sys
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Make regions from bam header')
    parser.add_argument('--chunk_size', action="store", dest='chunk_size')
    parser.add_argument('--region', action="store", dest='region')
    parser.add_argument('--bam_header_path', action="store", dest='bam_header_path')
    
    args = parser.parse_args()
    chunk_size = args.chunk_size
    region = args.region
    header_path = args.bam_header_path

    if chunk_size is None and region is not None:
        print("Error: chunk_size is a required input if region is inputted")
        print("Usage: python3 chunky_regions.py --chunk_size INT --region [CHR:START-END] --bam_header_path [PATH]")
        print("or")
        sys.exit("Usage: bam_header stdin | python3 chunky_regions.py --chunk_size INT --region [CHR:START-END]")
    
    if chunk_size is not None:
        chunk_size = int(chunk_size)

    if header_path is None:
        # print("Getting the header from stdin as it was not provided as an input")
        if sys.stdin.isatty():
            print("Error: No bam header supplied to stdin")     
            print("Usage: python3 chunky_regions.py --chunk_size INT --region [CHR:START-END] --bam_header_path [PATH]")
            print("or")
            sys.exit("Usage: bam_header stdin | python3 chunky_regions.py --chunk_size INT --region [CHR[:START-END]]")
        else:
            bam_header = sys.stdin

    else:
        bam_header = open(header_path)                      

    # print(bam_header) 
    bam_header_contents = [line.rstrip('\n') for line in bam_header]
    # print(bam_header_contents)
    # print([line[0:3] for line in bam_header_contents])
    bam_header_contents = [line for line in bam_header_contents  if line[0:3] == "@SQ" ]
    # print(x)
    # for line in bam_header_contents:
        # print(line)
    bam_header_contents = [line.split("\t") for line in bam_header_contents]
    # print(bam_header_contents)
    bam_chrom = [line[1] for line in bam_header_contents]
    bam_chrom = [line.split("SN:")[1] for line in bam_chrom]
    length = [line[2] for line in bam_header_contents]
    length = [line.split("LN:")[1] for line in length]
    # print(bam_chrom)
    # print(length)

    if region is None:
        # print("No region supplied - will create chunks over entrire bam header")
        sys.stderr.write("No region supplied - will create chunks over entrire bam header \n")
        if chunk_size is not None:

            for chrom_name in bam_chrom:
                # print(chrom_name)
                chrom_idx = bam_chrom.index(chrom_name)
                chrom_length = int(length[chrom_idx])
                chunk_start = 0

                while chunk_start < chrom_length:
                    chunk_end = chunk_start + chunk_size
                    if chunk_end > chrom_length:
                        chunk_end = chrom_length
                    print(chrom_name + ":" + str(chunk_start) + "-" + str(chunk_end))
                    chunk_start = chunk_end
        if chunk_size is None:
            sys.stderr.write("No chunk_size supplied - will just output chromosomes \n")
            for chrom in bam_chrom:
                print(chrom)
            


    if region is not None:
        region = region.split(":")
        # print(region)
        if len(region) == 1:
            # split chromosome into equal size chunks
            chrom_name = region[0]
            chrom_idx = bam_chrom.index(chrom_name)
            chrom_length = int(length[chrom_idx])
            # print(chrom_length)
            chunk_start = 0

            while chunk_start < chrom_length:
                # start = region_start
                chunk_end = chunk_start + chunk_size
                if chunk_end > chrom_length:
                    chunk_end = chrom_length
                print(chrom_name + ":" + str(chunk_start) + "-" + str(chunk_end))
                chunk_start = chunk_end
        else:
            # split chr:start-end into equal size chunks
            chrom_name = region[0]
            start_end = region[1]
            start_end = start_end.split("-")
            if len(start_end) == 2:
                region_start = int(start_end[0])
                region_end = int(start_end[1])

                chrom_idx = bam_chrom.index(chrom_name)
                chrom_length = int(length[chrom_idx])
                # print(chrom_length)
                chunk_start = region_start
                chrom_length = min(chrom_length, region_end)
                while chunk_start < chrom_length:
                    # start = region_start
                    chunk_end = chunk_start + chunk_size
                    if chunk_end > chrom_length:
                        chunk_end = chrom_length
                    print(chrom_name + ":" + str(chunk_start) + "-" + str(chunk_end))
                    chunk_start = chunk_end

            else:
                sys.exit("Error: region must be CHR:START-END or just CHR")


    