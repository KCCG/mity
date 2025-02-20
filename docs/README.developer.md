<!-- omit in toc -->
# Developer Guide

- [Local development environment](#local-development-environment)
  - [Venv](#venv)
  - [Pip compile locally](#pip-compile-locally)
  - [Local testing](#local-testing)
- [Development Branch](#development-branch)
  - [TestPyPI Repo](#testpypi-repo)
  - [Test DockerHub Repo](#test-dockerhub-repo)
- [Pysam](#pysam)
  - [VariantFile](#variantfile)
  - [VariantHeader](#variantheader)
  - [VariantRecord](#variantrecord)
  - [VariantRecordSamples](#variantrecordsamples)
- [MITY call](#mity-call)
  - [Freebayes](#freebayes)
- [MITY normalise](#mity-normalise)
  - [Bcftools norm](#bcftools-norm)
- [MITY report](#mity-report)
  - [vcfanno](#vcfanno)
  - [Vep annotations](#vep-annotations)
    - [Vep labels](#vep-labels)
    - [Vep values](#vep-values)
  - [Pandas and Excel](#pandas-and-excel)
- [Updating mitomap annotation sources](#updating-mitomap-annotation-sources)
  - [Conversion scripts](#conversion-scripts)
  - [MT and chrM](#mt-and-chrm)
  - [Updating vcfanno config](#updating-vcfanno-config)

## Local development environment

### Venv

### Pip compile locally

### Local testing

## Development Branch

Merging into the development branch triggers an Azure pipeline which performs the following:

- compiles the python package and uploads it to `TestPyPI`.
- creates and publishes a docker image.

Both repos will have a tag in the form (where `r` is the revision number):

```txt
YYYY.MMDD.r
```

### TestPyPI Repo

<https://test.pypi.org/project/mitywgs-test/>

### Test DockerHub Repo

The dockerhub test repo is currently private.

## Pysam

A large part of mity uses pysam, which has quite poor documentation. So here's a guide for mity's use cases. For more details, I would highly recommend going straight to the source code and trying things out to test behaviour.

<https://github.com/pysam-developers/pysam/blob/cdc0ed12fbe2d7633b8fa47534ab2c2547f66b84/pysam/libcbcf.pyx>

### VariantFile

```python
import pysam

variant_file_obj = pysam.VariantFile(vcf_file_name)
```

### VariantHeader

```python
header_obj = variant_file_obj.header
```

Adding headers:

```txt
header_obj.filters.add(params)
header_obj.info.add(params)
header_obj.formats.add(params)
```

### VariantRecord

To loop through variants:

```python
for variant in variant_file_obj.fetch():
  ## do something with variant
```

To access fields:

```python
variant.chrom
variant.pos
variant.ref
variant.alts  ## note that this is a tuple of all the alts
variant.qual
variant.filter
```

Note that `variant.filter` is a dictionary, so it should have both keys and values, but filter names don't have a value, so the keys are the names and the values are just blank.

```python
variant.filter.keys()   ## this is a tuple/list
```

To get info columns, we can use `variant.info` which is a dictionary.

```python
variant.info["info column name"]
```

Note that the dictionaries here support all the the dictionary syntax, e.g.

```python
for info_field_name, value in variant.info.items():
```

### VariantRecordSamples

To loop through samples:

```python
for sample in variant.samples.values():
  ## do something with sample
```

To get the sample name:

```python
sample.name
```

We can access sample fields like a dictionary:

```python
sample["field"]
```

## MITY call

### Freebayes

<https://github.com/freebayes/freebayes>

## MITY normalise

### Bcftools norm

Documentation:
<https://samtools.github.io/bcftools/bcftools.html#norm>

Bcftools command:

```python
pysam.bcftools.norm("-f", self.reference_fasta, "-m-both", self.vcf)
```

## MITY report

### vcfanno

More information can be found here: [vcfanno](https://github.com/brentp/vcfanno)

Since vcfanno tends to have long, somewhat verbose warnings, we capture the `stdout` and it's only displayed in debug mode.

```python
    res = subprocess.run(
        vcfanno_cmd,
        shell=True,
        check=False,
        ## capture the output of the vcfanno command since it tends to produce
        ## long warning messages, output shown in --debug mode
        capture_output=True,
        text=True,
    )

    logger.debug("vcfanno output:")
    logger.debug(res.stdout)
```

### Vep annotations

#### Vep labels

Example header INFO line:

```text
##INFO=<ID=CSQ,Number=.,Type=String,Description=
"Consequence annotations from Ensembl VEP. Format:
Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...
```

Returns a list of keys, i.e.

```python
["Allele", "Consequence", ...]
```

Change this line if the description text or format changes.

```python
description = description.replace(
    "Consequence annotations from Ensembl VEP. Format:", ""
)
```

#### Vep values

In function `get_vep_values(self, variant)`.

```text
Takes a string from VEP consequences/impacts in the form:
    impact value | impact value | ... |, |||, |||

And peforms the following:
    - removes annotation "line" if the consequence has "stream" in it,
      e.g. "upstream"
    - concatenates remaining line fields with ";"

Types:
    variant.info["CSQ"]: Tuple of form (a|b|c|..., a|b|c|..., a|b|c|...)

Example:
    vep_keys = [ IMPACT, Consequence, field_1, field_2, field_3 ]

    list(variant.info["CSQ"]) = [
        HIGH        | something_else    | a | b | c,
        LOW         | upstream_variant  | d | e | f,
        MODIFIER    | something_else    | g | h | i
    ]

    vep_dict = {
        HIGHEST IMPACT VEP  : HIGH,
        CONSEQUENCE VEP     : something_else;something_else,
        FIELD 1 VEP         : a;g,
        FIELD 2 VEP         : b;h,
        FIELD 3 VEP         : c;i
    }
```

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

## Updating mitomap annotation sources

### Conversion scripts

Many of the tables and resources from Mitomap are not in vcf format, which we need to annotate using vcfanno.

### MT and chrM

To support both notations, we create two sets of the same annotation file: one with `MT` and the other with `chrM`. Aside from that, these annotation files are identical. To generate the `chrM` versions, you can simply find and replace `MT` and `chrM`.

### Updating vcfanno config

Vcfanno config files live in `mitylib/config/` as `vcfanno-config-chrm.toml` and `vcfanno-config-mt.toml` for `chrM` and `MT` respectively.

Refer to: [brentp/vcfanno](https://github.com/brentp/vcfanno) for more detailed information about vcfanno. But as a primer, this is the basic layout of an annotation:

```toml
[[annotation]]
file="path/to/file.vcf"
fields=["a", "b"] # field names in the vcf file
names=["renamed_a", "renamed_b"] # renamed field names
ops=["self", "self"] # ops=operations, in this case self = don't do anything
```
