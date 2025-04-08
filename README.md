![mity logo](https://github.com/KCCG/mity/blob/master/logos/mity-logo-red-white.png "mity")

# mity

`mity` is a bioinformatic analysis pipeline designed to call mitochondrial SNV and INDEL variants from Whole Genome Sequencing (WGS) data. `mity` can:

* identify very low-heteroplasmy variants, even <1% heteroplasmy when there is sufficient read-depth (eg >1000x)
* filter out common artefacts that arise from high-depth sequencing
* easily integrate with existing nuclear DNA analysis pipelines (mity merge)
* provide an annotated report, designed for clinicians and researchers to interrogate

# Usage

```bash
mity -h
```

More detailed usage can be found [docs/commands.md](docs/commands.md)

# Dependencies

* python3 (tested on 3.10, 3.12)
* freebayes >= 1.2.0
* bgzip + tabix
* gsort (<https://github.com/brentp/gsort>)
* pyvcf
* xlsxwriter
* pandas

# Installation

Installation instructions via Docker, pip, or manually are available in [INSTALL.md](docs/INSTALL.md)

# Example Usage

This is an example of calling variants in the Ashkenazim Trio.

## mity call

First run `mity call` on three MT BAMs provided in [docs/test_files.md](docs/test_files.md). CRAM files are supported.

We recommend always using `--normalise`, or `mity report` won't work:

```bash
mity call \
--prefix ashkenazim \
--output-dir test_out \
--region MT:1-500 \
--normalise \
test_in/HG002.hs37d5.2x250.small.MT.RG.bam \
test_in/HG003.hs37d5.2x250.small.MT.RG.bam \
test_in/HG004.hs37d5.2x250.small.MT.RG.bam 
```

This will create `test_out/normalised/ashkenazim.mity.vcf.gz` (and tbi file).

or, if using Docker:

```bash
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity call \
--prefix ashkenazim \
--output-dir test_out \
--region MT:1-500 \
--normalise \
input/HG002.hs37d5.2x250.small.MT.RG.bam \
input/HG003.hs37d5.2x250.small.MT.RG.bam \
input/HG004.hs37d5.2x250.small.MT.RG.bam 
```

## mity normalise

High-depth sequencing and sensitive variant calling can create many variants with more than 2 alleles, and in some
cases, joins two nearby variants separated by shared `REF` sequence into a multi-nucleotide polymorphism
as discussed in the manuscript. Here, variant normalisation relates to decomposing the multi-allelic variants and
where possible, splitting multi-nucleotide polymorphisms into their cognate smaller variants. At the time of writing,
all variant decomposition tools we used failed to propagate the metadata in a multi-allelic variant to the split
variants which caused problems when reporting the quality scores associated with each variant.
  
Technically you can run `mity call` and `mity normalise` separately, but since `mity report` requires a normalised
vcf file, we recommend running `mity call --normalise`.

```bash
mity normalise \
--prefix ashkenazim \
--output-dir output \
output/ashkenazim.mity.call.vcf.gz
```

## mity report

We can create a `mity report` on the normalised VCF:

```bash
mity report \
--prefix ashkenazim \
--min_vaf 0.01 \
--output-dir output \
output/ashkenazim.mity.normalise.vcf.gz
```

This will create: `output/ashkenazim.mity.report.xlsx`.

For more information about creating custom report configurations, see [custom_report_configs](docs/custom_report_configs.md).

## mity merge

You can merge a nuclear vcf.gz file and a mity.vcf.gz file thereby replacing the MT calls from the nuclear VCF (
presumably from a caller like HaplotypeCaller which is not able to sensitively call mitochondrial variants) with
the calls from `mity`.

```bash
mity merge \
--prefix ashkenazim \
--mity_vcf output/ashkenazim.mity.vcf.gz \
--nuclear_vcf todo-create-example-nuclear.vcf.gz \
--output-dir output
```

## mity runall

To run `call`, `normalise` and `report` all in one go, you can use the `mity runall` command. This command supports all the options from `call`, `normalise` and `report`.

## Other usage information

All commands have the options for debugging and analysis:

```bash
-k, --keep            Keep all intermediate files
-d, --debug           Enter debug mode
```

All commands have the options for consistent naming:

```bash
--prefix PREFIX       Output files will be named with PREFIX
--output-dir OUTPUT_DIR
                        Output files will be saved in OUTPUT_DIR. Default: '.'
```

All commands generate output files in a structured manner:

| Mity Command   | Files Generated                                   |
| -------------- | ------------------------------------------------- |
| mity call      | prefix.mity.call.vcf.gz                           |
| mity normalise | prefix.mity.normalise.vcf.gz                      |
| mity report    | prefix.mity.annotated.vcf (with --keep)<br>prefix.mity.report.xlsx |
| mity merge     | prefix.mity.merge.vcf.gz                          |

# Recommendations for interpreting the report

For more information about report columns and annotation sources, see [mity_report_documentation](docs/mity_report_documentation.md).

Assuming that you are looking for a pathogenic variant underlying a patient with a rare genetic disorder potentially
caused by a Mitochondrial mutation, then we recommend the following strategy:

1. tier 1 or 2 variants included in the 'commercial_panels' column
2. tier 1 or 2 variants that match the clinical presentation and the phenotype in 'disease_mitomap', preferably
those that are annotated with Confirmed evidence in the 'status_mitomap' column
3. exclude common variants: anything linked to 'phylotree_haplotype', high 'phylotree_haplotype', high
'MGRB_frequency', high 'GenBank_frequency_mitomap'.
4. consider any remaining tier 1 or 2 variants that may have a predicted impact on tRNA
5. consider any remaining variants with high numbers of 'variant_references_mitomap'
6. if you have analysed multiple family members, consider variants who's level of 'variant_heteroplasmy' match the
disease burden
7. tier 3 variants have low numbers of supporting reads, and should be considered with caution. However we have observed
numerous tier 3 variants, especially in WGS from blood, that match the pathogenic allele known to be at much higher
heteroplasmy in the affected tissue (this phenomenon is well established in the literature). Thus, if there are any
tier 3 variants identified that match the patient's clinical presentation, then we recommend considering these
as candidate variants and validating using an orthogonal clinically validated assay, preferably on the disease
affected tissue.

# Reference genomes

## Human

`mity` natively supports the analysis of the revised Cambridge Reference Sequence (rCRS, RefSeq ID NC_012920.1). The
rCRS used in most human reference genomes from NCBI (GRCh37, hs37d5, GRCh38) and hg38 from UCSC, where it is either
named `chrM`, or `MT`. The main exception in common use is the `hg19` reference genome from UCSC, which used a different
sequence (RefSeq NC_001807) which differs in length by 2bp, and sharing 99% sequence homology (16530/16572 identities)
and 4 gaps. For now, `mity call` supports the hg19 reference, but `mity report` will not annotate variants properly, so
you should not use this part of the pipeline. We strongly recommend that for mitochondrial analysis, to use a reference
genome that uses the rCRS sequence.

> * the mitochondrial genome: since the release of the UCSC hg19
> assembly, the Homo sapiens mitochondrion sequence (represented as "chrM" in the
> Genome Browser) has been replaced in GenBank with the record NC_012920, the
> revised Cambridge Reference Sequence (rCRS).  We have not replaced the original
> sequence, NC_001807, as chrM in the hg19 Genome Browser.  However, files in the
> subdirectory p13.plusMT include NC_012920 as "chrMT", in addition to the original
> "chrM".

| Reference   | contig name | RefSeq ID   | length   | rCRS |
| ----------- | ----------- | ----------- | ---------| ---- |
| GRCh37      | chrM        | NC_012920.1 | 16569 bp | rCRS |
| hs37d5      | MT          | NC_012920.1 | 16569 bp | rCRS |
| hg19 (UCSC) | chrM        | NC_001807.4 | 16571 bp | no   |
| GRCh38      | chrM        | NC_012920.1 | 16569 bp | rCRS |

## Mouse

`mity` `call` and `normalise` support the analysis of the mouse genome (`mity call --reference mm10 ...`). `mity report`
currently only supports variant annotation to the human rCRS sequence.

## Custom reference fasta and genome files

`mity call`, `normalise` and `runall` support custom reference fasta and genome files with the options `--custom-reference-genome`
 and `--custom-reference-fasta`.

# Commonly asked Questions

## Base quality score recalibration (BQSR)

Most of the development of `mity` was tested on BAM files that had undergone GATK's BQSR method, which improves the
base qualities of each read.
In our experience, this reduced the quality score of most bases by ~10 points, indicating that the base qualities
straight out of the sequencer are generally inflated. As the GATK best practices guide no longer recommends BQSR, it's
reasonable to ask whether `mity` can be run on BAM files straight out of the aligner.
`mity` has a custom QUAL score, which depends on the base qualities of only the reads that support the alternative
allele.  
For tier 1 or 2 variants, there will be so many supporting reads, that any miscalibration of base quality scores will
have no material effect. Tier 3 variants with very few supporting reads may be impacted, where a variant with only 3 or
4 supporting reads may end up having a stronger mity QUAL score than after BQSR. The comment above regarding how you
should interpret and validate tier 3 variant still holds.
We would appreciate any feedback you may have on this.

## CRAM support

CRAM support was added to `mity call` in v0.4.0.

# Acknowledgements

We would like to thank:

* The Kinghorn Centre for Clinical Genomics and collaborators, who helped with feedback for running `mity`.
* The Genome in a Bottle consortium for providing the test data used here
* Eric Talevich who's CNVkit helped us structure `mity` as a package
* Erik Garrison for developing `FreeBayes` and his early feedback in optimising `FreeBayes` for sensitive variant detection.
* Brent Pederson for developing `gsort`
