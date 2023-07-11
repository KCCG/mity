# Test Scripts

```
docker run -w "$PWD" -v "$PWD":"$PWD" mity-test-0.4.7 call \
--prefix ashkenazim \
--out-folder-path tests/ashkenazim-trio/output \
--region MT:1-500 \
--normalise \
tests/ashkenazim-trio/input/HG002.hs37d5.2x250.small.MT.RG.bam \
tests/ashkenazim-trio/input/HG003.hs37d5.2x250.small.MT.RG.bam \
tests/ashkenazim-trio/input/HG004.hs37d5.2x250.small.MT.RG.bam 
```