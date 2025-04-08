# change log

## 0.4.0

- Includes x38 support, CRAM, mitomap scores and pathogenic variant annotations from MitoTip.

## 1.0.0

- Refactored all modules of mity to support newer versions of python and other dependencies.
- Updated mity to use pysam for all vcf / bam / cram file manipulation.
- Updated mity report to use vcfanno for annotating variants.
- Updated mity merge to use pysam concat and pysam isec for merging the files.
- General performance enhancements across all mity modules.

### 1.1.0

- Fixed an issue where the `report-config.yaml` and `vcfanno-config.toml` files where not being found by mity.
- Added options to specify custom report and vcfanno config files.
- Added `--contig` as an option to `mity report` to specify the contig.
- Added annotations for `MT` and `chrM`. `mity` will default to using `MT` if the `--contig` option is not provided.
- Added vcfanno as a dependency to the `Dockerfile`.
- Added a `runall` option to run `mity call`, `mity normalise`, and `mity report` sequentially. This option supports all options from `call` and `report`.

### 1.2.0

- Added new option for `mity report` to additionally output the annotated vcf file as a xlsx file.
- Added and improved documentation including: [change_log](./change_log.md), [custom_report_configs](./custom_report_configs.md), [mity_report_documentation](./mity_report_documentation.md), [future_improvements](./future_improvements.md).
- Fixed but in mity report when the sheet name is longer that 31 characters.

## 2.0.0

- Updated and added new mitomap sources (see [mity_report_documentation](./mity_report_documentation.md) for more detailed updates). Note that these changes are breaking and affect the behaviour of `mity report`.
- Added new option for `mity report` to additionally output the annotated vcf file as a xlsx file.
- Fixed `mity report` output file name from `prefix.annotated_variants.xlsx -> prefix.mity.report.xlsx`.
- Added and improved documentation including: [change_log](./change_log.md), [custom_report_configs](./custom_report_configs.md), [mity_report_documentation](./mity_report_documentation.md), [future_improvements](./future_improvements.md).
- Added `ruff.toml` for formatting.
- Fixed but in mity report when the sheet name is longer that 31 characters.
