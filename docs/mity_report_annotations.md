# mity report annotations

If you are using the pip installation method, then you can use custom annotations.
Default annotations can be found in the `mitylib` directory.

## Adding new annotations

1. Place annotations file (and index file) in `mitylib/annot`.
2. Add a new vcfanno entry for the file in `vcfanno-config.toml`
3. Add new fields in `annotations` and `excel_headers` in `report-config.yaml`.

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

Note that `your_annotation_field_name` needs to match the info field name exactly.
