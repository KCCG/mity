<!-- omit in toc -->
# Developer Guide

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
- [Future Improvements](#future-improvements)
  - [General Options](#general-options)
    - [MITY Run](#mity-run)
    - [Run on all files in directory](#run-on-all-files-in-directory)
  - [Multithreading](#multithreading)
    - [Bcftools](#bcftools)
    - [MITY operation threading](#mity-operation-threading)
      - [Python Example](#python-example)
      - [Possible MITY Implementation](#possible-mity-implementation)
  - [Cache (with lru\_cache)](#cache-with-lru_cache)
    - [Python lru\_cache Example (fibonacci)](#python-lru_cache-example-fibonacci)
    - [Possible MITY Implementation](#possible-mity-implementation-1)
  - [MITY CLI and Commands](#mity-cli-and-commands)
  - [MITY Report](#mity-report-1)
    - [Specify default annotation sources](#specify-default-annotation-sources)
    - [Custom annotation sources](#custom-annotation-sources)
    - [Move annotation sources to s3](#move-annotation-sources-to-s3)


# Development Branch

Merging into the development branch triggers an Azure pipeline which performs the following:
- compiles the python package and uploads it to `TestPyPI`.
- creates and publishes a docker image.

Both repos will have a tag in the form (where `r` is the revision number):
```
YYYY.MMDD.r
```

## TestPyPI Repo

<https://test.pypi.org/project/mitywgs-test/>

## Test DockerHub Repo

The dockerhub test repo is currently private.

# Pysam

A large part of mity uses pysam, which has quite poor documentation. So here's a guide for mity's use cases. For more details, I would highly recommend going straight to the source code and trying things out to test behaviour.

<https://github.com/pysam-developers/pysam/blob/cdc0ed12fbe2d7633b8fa47534ab2c2547f66b84/pysam/libcbcf.pyx>

## VariantFile

```python
import pysam

variant_file_obj = pysam.VariantFile(vcf_file_name)
```

## VariantHeader

```python
header_obj = variant_file_obj.header
```

Adding headers:
```
header_obj.filters.add(params)
header_obj.info.add(params)
header_obj.formats.add(params)
```

## VariantRecord

To loop through variants:
```python
for variant in variant_file_obj.fetch():
  # do something with variant
```

To access fields:
```python
variant.chrom
variant.pos
variant.ref
variant.alts  # note that this is a tuple of all the alts
variant.qual
variant.filter
```

Note that `variant.filter` is a dictionary, so it should have both keys and values, but filter names don't have a value, so the keys are the names and the values are just blank.
```python
variant.filter.keys()   # this is a tuple/list
```

To get info columns, we can use `variant.info` which is a dictionary.
```python
variant.info["info column name"]
```

Note that the dictionaries here support all the the dictionary syntax, e.g.
```python
for info_field_name, value in variant.info.items():
```

## VariantRecordSamples

To loop through samples:
```python
for sample in variant.samples.values():
  # do something with sample
```

To get the sample name:
```python
sample.name
```

We can access sample fields like a dictionary:
```python
sample["field"]
```


# MITY call

## Freebayes

<https://github.com/freebayes/freebayes>

# MITY normalise

## Bcftools norm

Documentation:
<https://samtools.github.io/bcftools/bcftools.html#norm>

Bcftools command:
```
pysam.bcftools.norm("-f", self.reference_fasta, "-m-both", self.vcf)
```


# MITY report

## vcfanno

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

## Vep annotations

### Vep labels

Example header INFO line:
```
##INFO=<ID=CSQ,Number=.,Type=String,Description=
"Consequence annotations from Ensembl VEP. Format:
Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature|...
```

Returns a list of keys, i.e. 
```
["Allele", "Consequence", ...]
```

Change this line if the description text or format changes.
```
description = description.replace(
    "Consequence annotations from Ensembl VEP. Format:", ""
)
```

### Vep values

In function `get_vep_values(self, variant)`.

```
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

## Pandas and Excel

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

# Future Improvements

## General Options

### MITY Run

Option to run mity call, normalise and report all in one go.

```python
from mitylib import call, normalise, report

class MityRun:
	def __init__(self, **all_params_needed_for_call_norm_report):
		# init variables here
		self.run()

	def run(self):
		call.Call(**params)
		normalise.Normalise(**params)
		report.Report(**params)
```

Changes needed:
- `commands.py` changes to include a new subparser `P_run`
- New `run.py` module
### Run on all files in directory

Currently the only way to specify a file is to put it as a command line argument, but we can add an option to run mity commands on all files in a directory.

commands.py
```python
P_xxxxx.add_argument(
    "--input-dir", # or "--directory" but we already have --output-dir so...
    action="store",
    type=str,
    required=False,
    help="Directory for input files",
    dest="input_dir",
)
```

each-module.py
```python
def run(self):
	...
	if self.input_dir is not None:
		files = [f for f in os.listdir(directory) if os.path.isfile(os.path.join(directory, f))]
		for file in files:
			run_single(file, params)
```

Changes needed:
- `commands.py` extra argument for call, normalise, report (see code above)
- `normalise.py` changes to support multiple files
- Changes to each module to check for files in a directory (see code above)
- Option in `report.py` to put each report in a separate excel file, or to keep them as separate sheets in the same excel file

## Multithreading

### Bcftools
Some bcftools commands also support multithreading, which can be enabled with arguments. The commands we use that support multithreading include:
- Normalise
- Concat

### MITY operation threading

Similarly, some mity operations that operate on variants independently (e.g. filtering in mity normalise) can be done concurrently.

#### Python Example

```python
import threading

def f(args):
	print("hello-world")

if __name__ == "__main__":
	t1 = threading.Thread(target=f, args=(args,))
	t2 = threading.Thread(target=f, args=(args,))

	t1.start() # starts thread 1
	t2.start() # starts thread 2
	
	t1.join() # stops "main" thread until t1 is complete
	t2.join() # same for t2
```

#### Possible MITY Implementation

```python
import threading

def filter_variants(variant):
	# filtering logic
	return filtered_variant

def run_filtering(vcf_obj):
	for variant in vcf_obj.records.fetch():
		threading.Thread(TODO)
```

## Cache (with lru_cache)

Some parts of mity that perform the same operation multiple times, can be sped up using python's `lrucache`.

### Python lru_cache Example (fibonacci)

```python
@lru_cache(maxsize=128)
def fibonacci(n):
	if n < 2:
		return n
	return fibonacci(n-1) + finbonacci(n-2)
```

### Possible MITY Implementation

```python
@lru_cache
def clean_string(self, s):
	"""
	Removes the following characters from a string:
		"
		'
		()

	Used to output text from annotation sources.
	"""
	if isinstance(s, tuple):
		s = s[0]

	s = str(s)
	s = s.replace('"', "").replace("(", "").replace(")", "")

	return s
```

- Often times there are the same repeated strings that need to be "cleaned", and this can be sped up by just reading from the previous stored output.
- There are definitely other methods/functions in MITY that have this same problem, and they can all be improved with this.
## MITY CLI and Commands

- Click instead of argparse
- Poetry instead of setuptools
## MITY Report

Most of the changes and improvements here come from making annotations (and potentially default annotations) more flexible and extendable.

### Specify default annotation sources

Example usage:
```bash
$ mity report --annotations anticodon mitomap --vcf prefix.mity.normalise.vcf.gz
```

This should only include the specified annotations in the report.

Changes needed:
- Separate labels chosen in the yaml file for annotations into each file, we can still run all the default annotations, but just pick the headers that we want to include in the final excel file based on what was chosen.
- Make annotations a bit more flexible overall

These changes should help the rest of the improvements too.

### Custom annotation sources

Allow users to add custom annotation sources to the mity report output using `vcfanno`.

Example usage:
```bash
$ mity report append-annotation \
	--vcfanno-config <vcfanno-config-file.toml> \
	--annotation <annotation.vcf>
	--fields <fields.yaml>
	
$ mity report replace-annotation \
	--vcfanno-config <vcfanno-config-file.toml> \
	--annotation <annotation.vcf>
	--fields <fields.yaml>
```

- `append-annotation` keeps the original mity annotations and adds new ones
- `replace-annoation` only uses the custom annotation

**`vcfanno-config-file.toml`**
```
[[annotation]]
file="annotation.vcf"
columns=[X,X]
names=[name1, name2]
ops=["self", "self"]
```

**`fields.yaml` (option 1)**
```
annotation.vcf:
- name1
- name2
```

- Use default behaviour for creating mity report labels
- Read the yaml file in as a list

**`fields.yaml` (option 2)**
```
annotation.vcf:
	name1: CUSTOM_NAME_1
	name2: CUSTOM_NAME_2
```

- Use custom excel header names
- Read the yaml file as a dictionary

### Move annotation sources to s3

In a similar way to supporting custom annotation sources, `vcfanno` supports annotation sources from links, so we can potentially keep all our annotation sources on s3.

Changes needed:
- Upload annotation sources to Netapp/s3
- Update `vcfanno-config.toml` (i.e. the default annotation source) to use links