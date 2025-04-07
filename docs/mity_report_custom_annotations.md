# mity report custom annotations

## Custom report config

If you want to configure which annotations appear in the final report or rename any column titles, you can specify a custom report config used to generate the report. The existing `report-config.yaml` lives in `mitylib/config/`.

### Structure of `report-config.yaml`

We would recommend copying the existing report config as a base and editing it from there. The sections you should change are: `excel_headers` and `annotations`. If you change the other sections, you may break the custom calculations that appear in report.

```yaml
excel_headers:
- SAMPLE
- HGVS
- GENE/LOCUS
...
```

```yaml
vcf_headers:
  annotations:
    ANTICODON: ANTICODON
    GENE: GENE
    GENE_BIOTYPE: GENE BIOTYPE
    PHYLOTREE_MUT: PHYLOTREE MUT
    PHYLOTREE_HAPLOTYPE: PHYLOTREE HAPLOTYPE
    ...
```

If you would like to remove `PHYLOTREE HAPLOTYPE` from the report, remove the line in **both** `annotations` and `excel_headers`. Similarly if you want to add a label, add it to both `annotations` and `excel_headers`.

Note that you should only alter/remove the following lines in `excel_headers`:

```yaml
- ANTICODON
- GENE
- GENE BIOTYPE
- PHYLOTREE MUT
- PHYLOTREE HAPLOTYPE
- MGRB FILTER
- MGRB AN
- MGRB AC
- MGRB FREQUENCY
- MITOMAP DISEASE AC
- MITOMAP DISEASE AF
- MITOMAP DISEASE AACHANGE
- MITOMAP DISEASE HOMOPLASMY
- MITOMAP DISEASE HETEROPLASMY
- MITOMAP DISEASE PUBMED IDS
- MITOMAP DISEASE DISEASE
- MITOMAP DISEASE DISEASE STATUS
- MITOMAP DISEASE HGFL
- MITOMAP CONFIRMED MUTATIONS LOCUS
- MITOMAP CONFIRMED MUTATIONS LOCUSTYPE
- MITOMAP CONFIRMED MUTATIONS ASSOCIATEDDISEASE
- MITOMAP CONFIRMED MUTATIONS ALLELE
- MITOMAP CONFIRMED MUTATIONS AMINOACIDCHANGE
- MITOMAP CONFIRMED MUTATIONS STATUSMITOMAPCLINGEN
- MITOMAP CONFIRMED MUTATIONS LASTUPDATE
- MITOMAP MUTATIONS CODING CONTROL LOCUS
- MITOMAP MUTATIONS CODING CONTROL ALLELE
- MITOMAP MUTATIONS CODING CONTROL DISEASE
- MITOMAP MUTATIONS CODING CONTROL NUCLEOTIDECHANGE
- MITOMAP MUTATIONS CODING CONTROL AMINOACIDCHANGE
- MITOMAP MUTATIONS CODING CONTROL PLASMY
- MITOMAP MUTATIONS CODING CONTROL STATUS
- MITOMAP MUTATIONS CODING CONTROL GB FREQ
- MITOMAP MUTATIONS CODING CONTROL GB SEQS
- MITOMAP MUTATIONS CODING CONTROL REFERENCES
- MITOMAP MUTATIONS RNA LOCUS
- MITOMAP MUTATIONS RNA DISEASE
- MITOMAP MUTATIONS RNA ALLELE
- MITOMAP MUTATIONS RNA RNA
- MITOMAP MUTATIONS RNA HOMOPLASMY
- MITOMAP MUTATIONS RNA HETEROPLASMY
- MITOMAP MUTATIONS RNA STATUS
- MITOMAP MUTATIONS RNA MITOTIP
- MITOMAP MUTATIONS RNA GB FREQ
- MITOMAP MUTATIONS RNA GB SEQS
- MITOMAP MUTATIONS RNA REFERENCES
- MITOMAP POLYMORPHISMS AC
- MITOMAP POLYMORPHISMS AF
- MITOMAP POLYMORPHISMS HGFL
- MITOMAP VARIANTS CODING LOCUS
- MITOMAP VARIANTS CODING NUCLEOTIDECHANGE
- MITOMAP VARIANTS CODING CODONNUMBER
- MITOMAP VARIANTS CODING CODONPOSITION
- MITOMAP VARIANTS CODING AMINOACIDCHANGE
- MITOMAP VARIANTS CODING GB FREQ
- MITOMAP VARIANTS CODING GB SEQS
- MITOMAP VARIANTS CODING CURATEDREFERENCES
- MITOMAP VARIANTS CONTROL LOCUS
- MITOMAP VARIANTS CONTROL NUCLEOTIDECHANGE
- MITOMAP VARIANTS CONTROL GB FREQ
- MITOMAP VARIANTS CONTROL GB SEQS
- MITOMAP VARIANTS CONTROL CURATEDREFERENCES
- MITOTIP SCORE
- MITOTIP PERCENTILE
- MITOTIP QUARTILE
- MITOTIP SCORE INTERPRETATION
- MITOMAP STATUS
- COUNT
- PERCENTAGE
```

## Custom annotation sources

Adding custom annotation sources is a little tricker as the annotation sources are baked into mity. We would recommend:

1. Copying the annotation sources from the mity repo (i.e. this repo) from either the `annot_mt` or `annot_chrm` folders respectively. Put them into a local folder on your machine, e.g. `my_annotation_sources/`.
2. Add your own annotation sources to the same folder (ideally they're also gzipped and tabix'ed).
3. Copy the `vcfanno-config.yaml` file from the mity repo and use that as a base. Add your own annotation sources.
4. Copy the `report-config.yaml` file from the mity repo and use that as a base. Add you own labels for your own annotation sources.
5. When you run `mity report`, make sure to set:
   1. `--custom-vcf-config`
   2. `--custom-report-config`
   3. `--vcfanno-base-path`

Ensure that `vcfanno-config.toml` and `report-config.yaml` are in your current directory.

```yaml
excel_headers:
- SAMPLE
- HGVS
- GENE/LOCUS
...
- YOUR NEW EXCEL HEADER
```

Note that the order of headers in this yaml file represents the order of columns in the final excel file.

```yaml
vcf_headers:
  annotations:
    GenBank_frequency_mitomap: GENBANK FREQUENCY MITOMAP
    MGRB_AC: MGRB AC
    MGRB_AN: MGRB AN
    ...
    your_annotation_field_name: YOUR NEW EXCEL HEADER
```

Note that `your_annotation_field_name` needs to match the info field name from the `vcfanno` `names` field exactly.

See the existing config files as a guide for how to structure annotations.
