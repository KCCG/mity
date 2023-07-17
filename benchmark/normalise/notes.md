### differences

- mity's normalise is not left aligned, it does correctly separate differences though

### notes

- check whether the next step (mity report??) needs the exact output from mity normalise (i.e. can it use the output of bcftools normalise) 

### commands

```
bcftools norm -o bcftools_normalise.vcf -f ../../mitylib/reference/hs37d5.MT.fa -m-both ashkenazim.mity.vcf.gz
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity normalise -o mity_normalise.vcf.gz ashkenazim.mity.vcf.gz

vt normalize dbsnp.vcf -r seq.fa -o dbsnp.normalized.vcf
vt decompose_blocksub gatk.vcf -o decomposed_blocksub.vcf 
```
