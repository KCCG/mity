# Test Commands

Make sure to follow instructions in [INSTALL.md](INSTALL.md) to get the latest version of `mity`.

## Ashkenazim Trio

### Docker

MITY call and normalise:
```bash
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity call \
--prefix ashkenazim \
--out-folder-path <insert-output-folder-path> \
--region MT:1-500 \
--normalise \
ashkenazim-trio/input/HG002.hs37d5.2x250.small.MT.RG.bam \
ashkenazim-trio/input/HG003.hs37d5.2x250.small.MT.RG.bam \
ashkenazim-trio/input/HG004.hs37d5.2x250.small.MT.RG.bam 
```

MITY report:
```bash
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity report \
--prefix ashkenazim \
--min_vaf 0.01 \
--out-folder-path <insert-output-folder-path> \
<insert-output-folder-path>/ashkenazim.mity.vcf.gz
```

### PIP

MITY call and normalise:
```bash
mity call \
--prefix ashkenazim \
--out-folder-path <insert-output-folder-path> \
--region MT:1-500 \
--normalise \
ashkenazim-trio/input/HG002.hs37d5.2x250.small.MT.RG.bam \
ashkenazim-trio/input/HG003.hs37d5.2x250.small.MT.RG.bam \
ashkenazim-trio/input/HG004.hs37d5.2x250.small.MT.RG.bam 
```

MITY report:
```bash
mity report \
--prefix ashkenazim \
--min_vaf 0.01 \
--out-folder-path <insert-output-folder-path> \
<insert-output-folder-path>/ashkenazim.mity.vcf.gz
```