# future improvements

## General Options

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
mity report --annotations anticodon mitomap --vcf prefix.mity.normalise.vcf.gz
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

```toml
[[annotation]]
file="annotation.vcf"
columns=[X,X]
names=[name1, name2]
ops=["self", "self"]
```

**`fields.yaml` (option 1)**

```txt
annotation.vcf:
- name1
- name2
```

- Use default behaviour for creating mity report labels
- Read the yaml file in as a list

**`fields.yaml` (option 2)**

```txt
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
