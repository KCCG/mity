<!-- omit in toc -->
# Developer Guide

- [Development Branch](#development-branch)
  - [TestPyPI Repo](#testpypi-repo)
  - [Test DockerHub Repo](#test-dockerhub-repo)
- [Pysam VariantFile](#pysam-variantfile)
- [MITY call](#mity-call)
- [MITY normalise](#mity-normalise)
- [MITY report](#mity-report)
  - [vcfanno](#vcfanno)
  - [Dictionaries](#dictionaries)
  - [Cohort count and frequency](#cohort-count-and-frequency)
  - [Pandas and Excel](#pandas-and-excel)


## Development Branch

### TestPyPI Repo

### Test DockerHub Repo

## Pysam VariantFile


## MITY call
## MITY normalise
## MITY report

Structure:
1. call `vcfanno` to add annotations to the `mity.normalise.vcf.gz` file
2. parse the annotated vcf and organise into table with relevant information
3. convert to pandas df and finally write to excel spreadsheet

### vcfanno

More information can be found here: [vcfanno](https://github.com/brentp/vcfanno)

Since vcfanno tends to have long, somewhat verbose warnings, we capture the `stdout` and it's only displayed in debug mode.
```python
    res = subprocess.run(
        vcfanno_cmd,
        shell=True,
        check=False,
        # capture the output of the vcfanno command since it tends to produce
        # long warning messages, output shown in --debug mode
        capture_output=True,
        text=True,
    )

    logger.debug("vcfanno output:")
    logger.debug(res.stdout)
```

### Dictionaries

There are two main dictionaries in mity report:
- `headers_dict` contains all the columns for the excel spreadsheet in order. These are separated for readability, but they also correspond to sections of code. e.g. All the columns in `"start"` are under the comment:
```python
# headers_dict -> start
row.append(sample.name)
row.append(make_hgvs(variant.pos, variant.ref, variant.alts[0]))

for name in vcf_column_names["start_info"]:
    if name in variant.info.keys():
        row.append(clean_string(variant.info[name]))
        ...
```
- `vcf_column_names` contains the column names in the annotated vcf that correspond (mostly) to the columns in `headers_dict`. The columns that are missing are because that particular column does not directly relate to a field in the vcf and needs to be calculated e.g. HGVS, or is accessed in a special way, e.g. the position column.

### Cohort count and frequency

Since cohort count and frequency depends on the whole variant being read first, we add a placeholder at first, and fill it in after reading through all the samples.

```python
# placeholder for cohort count and cohort frequency
row.append(0)
row.append(0)
```

TODO: add other bit of code for counting and freq

### Pandas and Excel

We use the table and headers created in `make_table` and `make_headers` respectively to create the pandas dataframe:
```python
excel_pandasdf = pandas.DataFrame(excel_table, columns=excel_headers)
```
Then simply write to an excel file:
```python
xlsx_name = os.path.join(out_folder_path, prefix + ".annotated_variants.xlsx")
  with pandas.ExcelWriter(xlsx_name, engine="xlsxwriter") as writer:
      excel_pandasdf.to_excel(writer, sheet_name="Variants", index=False)
```