# Test Files

Test files are hosted on Netapp as zip files, and can be accessed via S3 commands.

To make repeated calls easier, you can place this alias in your `~/.bashrc` or `~/.zshrc`.

```bash
alias netapp="aws s3api --endpoint-url https://s3.ccia.org.au/ "
```

## Finding test files

To get a list of available test files, run this command:

```bash
netapp list-objects --bucket mity --prefix tests/ --query "Contents[].Key"
```

Or without the alias:

```bash
aws s3api --endpoint-url https://s3.ccia.org.au/ list-objects --bucket mity --prefix tests/ --query "Contents[].Key"
```

The output will be a list of keys that look like this:

```bash
[
    "tests/ashkenazim-trio.zip",
    "tests/bad-header.zip",
    "tests/cram-support.zip"
]
```

## Downloading test files

To download a test zip file, run this command:

```bash
netapp get-object --bucket mity --key <insert-key-here>
```

Or without the alias:

```bash
aws s3api --endpoint-url https://s3.ccia.org.au/ get-object --bucket mity --key <insert-key-here>
```

## Zip file structure

Each test zip file will have an `input` and `exp` folder. For example:

```bash
ashkenazim-trio
├── exp
│   ├── ashkenazim.mity.vcf.gz
│   └── ashkenazim.mity.vcf.gz.tbi
└── input
    ├── HG002.hs37d5.2x250.small.MT.RG.bam
    ├── HG002.hs37d5.2x250.small.MT.RG.bam.bai
    ├── HG003.hs37d5.2x250.small.MT.RG.bam
    ├── HG003.hs37d5.2x250.small.MT.RG.bam.bai
    ├── HG004.hs37d5.2x250.small.MT.RG.bam
    └── HG004.hs37d5.2x250.small.MT.RG.bam.bai
```

Specific commands to run each test can be found in [test-commands.md](test-commands.md)
