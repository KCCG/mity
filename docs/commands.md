# MITY Commands

If you're using `docker`, replace `mity` with:

```bash
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity
```

```bash
$ mity -h
usage: mity [-h] {call,normalise,report,merge,version} ...

Mity: a sensitive variant analysis pipeline optimised for WGS data

positional arguments:
  {call,normalise,report,merge,version}
                        mity sub-commands (use with -h for more info)
    call                Call mitochondrial variants
    normalise           Normalise & FILTER mitochondrial variants
    report              Generate mity report
    merge               Merging mity VCF with nuclear VCF
    version             Display this program's version.

options:
  -h, --help            show this help message and exit

Usage: See the online manual for details: http://github.com/KCCG/mity
Authors: Clare Puttick, Mark Cowley
License: MIT
Version: 0.4.45
```

## Call

```bash
usage: mity call [-h] [-d] [--reference {hs37d5,hg19,hg38,mm10}] [--prefix PREFIX] [--min-mapping-quality MIN_MQ] [--min-base-quality MIN_BQ] [--min-alternate-fraction MIN_AF] [--min-alternate-count MIN_AC] [--p P] [--normalise]
                 [--output-dir OUTPUT_DIR] [--region REGION] [--bam-file-list] [-k]
                 files [files ...]

positional arguments:
  files                 BAM / CRAM files to run the analysis on. If --bam-file-list is included, this argument is the file containing the list of bam/cram files.

options:
  -h, --help            show this help message and exit
  -d, --debug           Enter debug mode
  --reference {hs37d5,hg19,hg38,mm10}
                        Reference genome version to use. Default: hs37d5
  --prefix PREFIX       Output files will be named with PREFIX
  --min-mapping-quality MIN_MQ
                        Exclude alignments from analysis if they have a mapping quality less than MIN_MAPPING_QUALITY. Default: 30
  --min-base-quality MIN_BQ
                        Exclude alleles from analysis if their supporting base quality is less than MIN_BASE_QUALITY. Default: 24
  --min-alternate-fraction MIN_AF
                        Require at least MIN_ALTERNATE_FRACTION observations supporting an alternate allele within a single individual in the in order to evaluate the position. Default: 0.01, range = [0,1]
  --min-alternate-count MIN_AC
                        Require at least MIN_ALTERNATE_COUNT observations supporting an alternate allele within a single individual in order to evaluate the position. Default: 4
  --p P                 Minimum noise level. This is used to calculate QUAL score. Default: 0.002, range = [0,1]
  --normalise           Run mity normalise the resulting VCF
  --output-dir OUTPUT_DIR
                        Output files will be saved in OUTPUT_DIR. Default: '.'
  --region REGION       Region of MT genome to call variants in. If unset will call variants in entire MT genome as specified in BAM header. Default: Entire MT genome.
  --bam-file-list       Treat the file as a text file of BAM files to be processed. The path to each file should be on one row per bam file.
  -k, --keep            Keep all intermediate files
```

## Normalise

```bash
usage: mity normalise [-h] [-d] [--output-dir OUTPUT_DIR] [--prefix PREFIX] [--allsamples] [-k] [--p P] [--reference {hs37d5,hg19,hg38,mm10}] vcf

positional arguments:
  vcf                   vcf.gz file from running mity

options:
  -h, --help            show this help message and exit
  -d, --debug           Enter debug mode
  --output-dir OUTPUT_DIR
                        Output files will be saved in OUTPUT_DIR. Default: '.'
  --prefix PREFIX       Output files will be named with PREFIX
  --allsamples          PASS in the filter requires all samples to pass instead of just one
  -k, --keep            Keep all intermediate files
  --p P                 Minimum noise level. This is used to calculate QUAL scoreDefault: 0.002, range = [0,1]
  --reference {hs37d5,hg19,hg38,mm10}
                        Reference genome version to use. default: hs37d5
```

## Report

```bash
usage: mity report [-h] [-d] [--prefix PREFIX] [--min_vaf MIN_VAF] [--output-dir OUTPUT_DIR] [-k] vcf [vcf ...]

positional arguments:
  vcf                   mity vcf files to create a report from

options:
  -h, --help            show this help message and exit
  -d, --debug           Enter debug mode
  --prefix PREFIX       Output files will be named with PREFIX
  --min_vaf MIN_VAF     A variant must have at least this VAF to be included in the report. Default: 0.
  --output-dir OUTPUT_DIR
                        Output files will be saved in OUTPUT_DIR. Default: '.'
  -k, --keep            Keep all intermediate files
```

## Merge

```bash
usage: mity merge [-h] --mity_vcf MITY_VCF --nuclear_vcf NUCLEAR_VCF [--output-dir OUTPUT_DIR] [--prefix PREFIX] [--reference {hs37d5,hg19,hg38,mm10}] [-d] [-k]

options:
  -h, --help            show this help message and exit
  --mity_vcf MITY_VCF   mity vcf file
  --nuclear_vcf NUCLEAR_VCF
                        nuclear vcf file
  --output-dir OUTPUT_DIR
                        Output files will be saved in OUTPUT_DIR. Default: '.'
  --prefix PREFIX       Output files will be named with PREFIX. The default is to use the nuclear vcf name
  --reference {hs37d5,hg19,hg38,mm10}
                        reference genome version to use. default: hs37d5
  -d, --debug           Enter debug mode
  -k, --keep            Keep all intermediate files
```
