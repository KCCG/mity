# MITY-17 benchmark mity normalise vs bcftools

# research

### mitomap


### mity paper

> mity includes three additional filters: a strand bias filter to exclude variants with >90% or <10% alternative reads from
one strand, a read depth filter set to <15× depth, and a region filter to exclude variants in the homopolymeric regions at
m.302-319, or m.3105-3109, where there is an ‘N’ at m.3107, in the GRCh37 version of the mitochondrial sequence

### homopolymeric regions

>  Different sequencing platforms generate different types of sequencing errors, which can cause incorrectly called variants. The most common source of sequencing errors across platforms is the determination of nucleotides in so-called homopolymeric regions. These are regions that include stretches of the same nucleotide (e.g. AAAAA or TTTTTTTT). As a result of the internal chemistry used on platforms such as 454 and Ion Torrent, the number of identical nucleotides in such regions is often not accurately reported. This causes variant-callers to identify within homopolymer regions, insertions and deletions not actually present in the sample. The Illumina platform has a similar problem in which one nucleotide is surrounded by other nucleotides of the same type (e.g. AAAAGAAAA). Such cases are sometimes misread, with the different base identified as being the same as the surrounding nucleotides. This can lead to incorrect SNV calls. For example, a region of AAAAGAAAA in the sample may appear as AAAAAAAAA in the read. This could lead to a variant allele, A, being called where the G appears in the reference, when in fact the sample itself did contain a G at that position.
> 
> The Probabilistic Variant Caller includes an internal filter to recognize and prevent variants being reported in homopolymeric regions. 

Source: <https://resources.qiagenbioinformatics.com/manuals/clcgenomicsworkbench/650/Filtering_variants_in_homopolymeric_regions.html>

# mity normalise

## filtering

### code explanation

- Lines in the vcf file that don't pass have their filter columns annotated with the things that didn't pass. e.g. 

```
#CHROM	POS	ID	REF	ALT	QUAL	FILTER
MT	309	.	CTCCCCCGCTTCT	CCTCCCCCCGCTTCT	3220.0	SBA;POS
```

- `samp_X_fil` is an array containing 0's where the number of zeroes represent the number of samples that have failed that particular check. Then for each check, only one sample has to pass (it doesn't have to be the same sample that passes all the tests) for it to be given a PASS in the filter.


Format names
```
GT:DP:AD:RO:QR:AO:QA:GL
```

Format values
```
GT:DP:AD:RO:QR:AO:QA:GL	0/0:17028:17024,2:17024:653369:2:66:0,-5119.41,-58792.2	0/0:20957:20955,1:20955:803275:1:27:0 -6305.99,-72286.1	1/1:15676:5,15666:5:165:15666:598605:-53854.7,-4702.6,0
```

```python
for samp in FORMAT:
    
    samp_FORMAT = samp.split(":")

    AO = int(samp_FORMAT[FORMAT_names.index("AO")])
    RO = int(samp_FORMAT[FORMAT_names.index("RO")])
    AQR = float(samp_FORMAT[FORMAT_names.index("AQR")])

    from .util import make_hgvs
    hgvs = make_hgvs(int(line[1]), line[3], line[4])

    if RO > min_DP:
        # logging.debug("Position {} pass min_DP filter, with RO {} > {}".format(hgvs, str(RO), str(min_DP)))
        if not SB_range[0] <= SBR <= SB_range[1]:
            logging.debug("Position {} failed SBR filter {}".format(hgvs, str(SBR)))
            samp_SBR_fil.append(0)

        if MQMR < min_MQMR:
            samp_MQMR_fil.append(0)

        if AQR < min_AQR:
            samp_AQR_fil.append(0)

    if AO > min_DP:
        # logging.debug("Position {} pass min_DP filter, with AO {} > {}".format(hgvs, str(AO), str(min_DP)))
        if not SB_range[0] <= SBA <= SB_range[1]:
            logging.debug("Position {} failed SBA filter {}".format(hgvs, str(round(SBA, 3))))
            samp_SBA_fil.append(0)

    if POS in blacklist:
        samp_POS_fil.append(0)
```

```python
# if they aren't all zero than at least one passed
# only one sample has to pass
FILTER = []
if len(samp_SBR_fil) < len(FORMAT):
    SBR_fil = 1
else:
    SBR_fil = 0
    FILTER.append("SBR")

if len(samp_SBA_fil) < len(FORMAT):
    SBA_fil = 1
else:
    SBA_fil = 0
    FILTER.append("SBA")

if len(samp_MQMR_fil) < len(FORMAT):
    MQMR_fil = 1
else:
    MQMR_fil = 0
    FILTER.append("MQMR")

if len(samp_AQR_fil) < len(FORMAT):
    AQR_fil = 1
else:
    AQR_fil = 0
    FILTER.append("AQR")

if len(samp_POS_fil) < len(FORMAT):
    POS_fil = 1
else:
    POS_fil = 0
    FILTER.append("POS")

# FILTER = "FAIL"
if sum([SBR_fil, SBA_fil, MQMR_fil, AQR_fil, POS_fil]) == 5:
    FILTER = "PASS"
else:
    FILTER = ";".join(FILTER)
```


## quality value



# bcftools

## common options

```
-R      --regions-file FILE
```
> Regions can be specified either on command line or in a VCF, BED, or tab-delimited file (the default). The columns of the tab-delimited file can contain either positions (two-column format: CHROM, POS) or intervals (three-column format: CHROM, BEG, END), but not both. Positions are 1-based and inclusive. The columns of the tab-delimited BED file are also CHROM, POS and END (trailing columns are ignored), but coordinates are 0-based, half-open. To indicate that a file be treated as BED rather than the 1-based tab-delimited file, the file must have the ".bed" or ".bed.gz" suffix (case-insensitive). Uncompressed files are stored in memory, while bgzip-compressed and tabix-indexed region files are streamed. Note that sequence names must match exactly, "chr20" is not the same as "20". Also note that chromosome ordering in FILE will be respected, the VCF will be processed in the order in which chromosomes first appear in FILE. However, within chromosomes, the VCF will always be processed in ascending genomic coordinate order no matter what order they appear in FILE. Note that overlapping regions in FILE can result in duplicated out of order positions in the output. This option requires indexed VCF/BCF files. Note that -R cannot be used in combination with -r.

## norm

```
bcftools norm -o <outfile> -f <path_to_reference_fasta_file> -m-both <path_to_vcf.gz_file>
```

## filter

### options
```
-s      --soft-filter STRING|+
```
> annotate FILTER column with STRING or, with +, a unique filter name generated by the program ("Filter%d").

### commands

Removes blacklisted position (single)
```
bcftools filter --exclude "POS=310" <file>
```

Removes blacklisted position (mutliple)
```shell
bcftools filter --exclude "(POS > 301 && POS <319) || (POS >3104 && POS <3109)"  bcftools_normalise_test.vcf
```

Soft filter to replace the filter column with a string if it's excluded (does not remove the line)
```shell
bcftools filter --exclude "(POS > 301 && POS <319) || (POS >3104 && POS <3109)" --soft-filter "POS"  bcftools_normalise_test.vcf
```

# replacing normalise with pysam

## mini pysam.VariantFile.filter documentation

### add filter thing

1. Copy header from previous VCF file
2. Add filter options to header
3. Add filter to variant

```python
new_header = input_vcf.header
new_header.filters.add("FILTER_VALUE_1", None, None, "INSERT DESCRIPTION")
new_header.filters.add("FILTER_VALUE_2", None, None, "INSERT DESCRIPTION")
new_header.filters.add("FILTER_VALUE_3", None, None, "INSERT DESCRIPTION")

output_vcf = pysam.VariantFile('filtered.vcf', 'w', header=new_header)

BLACKLIST = list(range(302, 319)) + list(range(3105,3109))

for variant in input_vcf.fetch():
    # Adding one filter value
    variant.filter.add("FILTER_VALUE_1")

    # FILTER Output
    # FILTER_VALUE_1

    # Adding multiple filter values
    variant.filter.add("FILTER_VALUE_1")
    variant.filter.add("FILTER_VALUE_2")
    variant.filter.add("FILTER_VALUE_3")

    # FILTER Output
    # FILTER_VALUE_1;FILTER_VALUE_2;FILTER_VALUE_3

    output_vcf.write(variant)
```

Note that if you don't add a filter to your header, you get the following error:
```shell
File "pysam/libcbcf.pyx", line 2290, in pysam.libcbcf.VariantRecordFilter.add
KeyError: 'Invalid filter: pass'
```