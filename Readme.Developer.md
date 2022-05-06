
# Run version 1

## Call variants in the mitochondria
python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5 --prefix 15F00007 SYD-40100730.dedup.realigned.recalibrated.chrMT.bam SYD-40100741.dedup.realigned.recalibrated.chrMT.bam

## merge mity vcf and hc vcf
python3 merge.py --mity 15F00007.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --hc 15F00007.hc.small.vcf.gz

## make mity report
python3 report.py 15F00007.dedup.realigned.recalibrated.chrMT.mity.vcf.gz

###
# Run version 2
##

## Call variants in the mitochondria
python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5  SYD-40100741.dedup.realigned.recalibrated.chrMT.bam

python3 call.py --reference genome.MT.fa --min-alternate-fraction 0.5 SYD-40100730.dedup.realigned.recalibrated.chrMT.bam 

## make mity report
python3 report.py SYD-40100730.dedup.realigned.recalibrated.chrMT.mity.vcf.gz SYD-40100741.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --prefix 15F00007

# MJC
<!-- ## create fasta file
echo -e "MT\t1\t16569" > human_g1k_v37_decoy.MT.bed
samtools faidx $HOME/data/gatk-resource-bundle/2.8/b37/human_g1k_v37_decoy.fasta -o human_g1k_v37_decoy.MT.fa -r human_g1k_v37_decoy.MT.bed -->

## local tests
    mkdir test
    dx download project-F5bzbyQ0Jzv62xBg0P0Q5ygq:/test/kccg-mity-call/inputs/* -o test/
    freebayes -f /Users/marcow/data/gatk-resource-bundle/2.8/b37/human_g1k_v37_decoy.fasta -b test/A1.dedup.realigned.recalibrated.chrMT.bam -r MT:1-1000 --min-mapping-quality 30 --min-base-quality 20 --min-alternate-fraction 0.5 --min-alternate-count 4 --ploidy 2 --vcf unnormalised.vcf.gz
    mity call --reference $B37D5 --prefix test1 test/A1.dedup.realigned.recalibrated.chrMT.bam
    
    mity normalise --vcf test1.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --outfile test1.dedup.realigned.recalibrated.chrMT.mity.norm.vcf.gz
    ...
    ValueError: 'AD' is not in list

    mity report --vcf test1.dedup.realigned.recalibrated.chrMT.mity.vcf.gz 
    ...
    ValueError: 'SBR' is not in list

    dx ls project-F5bzbyQ0Jzv62xBg0P0Q5ygq:/test/kccg-mity-merge/inputs --brief | parallel dx download
    mity merge --mity_vcf test/15F00004.mity.vcf.gz --nuclear_vcf test/15F00004.hc.vqsr.vcf.gz --prefix xxx

### errors due to freebayes incompatibility?
local freebayes: version:  v1.0.2-dirty
kccg-freebayes: MD5 (resources/usr/bin/freebayes) = 88504cd29b834989f471cc119b748463
kccg-mity:      MD5 (resources/usr/bin/freebayes) = 88504cd29b834989f471cc119b748463
based on kccg-freebayes/Readme.Developer.md: this could be v1.0.2-33-gdbb6160 or maybe even 0.9.9

#### install freebayes-1.2.0
* This works great

    mkdir freebayes-1.2.0
    cd freebayes-1.2.0
    git clone --recursive --branch v1.2.0 git://github.com/ekg/freebayes.git
    make -j6
    export PATH=freebayes-1.2.0/freebayes/bin:$PATH
    which freebayes
    #freebayes-1.2.0/freebayes/bin/freebayes
    
    mity call --reference $B37D5 --prefix test1.2.0 test/A1.dedup.realigned.recalibrated.chrMT.bam --normalise
    mity call --reference $B37D5 --prefix test1.2.0 test/A1.dedup.realigned.recalibrated.chrMT.bam
    mity normalise --vcf test1.2.0.dedup.realigned.recalibrated.chrMT.mity.vcf.gz --outfile test1.2.0.dedup.realigned.recalibrated.chrMT.mity.norm.vcf.gz
    mity report --vcf test1.2.0.dedup.realigned.recalibrated.chrMT.mity.norm.vcf.gz
    mity merge --nuclear_vcf test/15F00004.hc.vqsr.vcf.gz --mity_vcf test/15F00004.mity.vcf.gz

# setup python package on pip
* https://packaging.python.org/tutorials/packaging-projects/
* I made setup.py based on this tutorial, and it seemed to work locally.
* use pip-compile from pip-tools to make requirements.txt: https://pypi.org/project/pip-tools/
* use build.sh to build sdist, bdist and then twine to upload to test pypi
* save twine credentials into keyring (pip install keyring)

    keyring set https://test.pypi.org/legacy/ drmjc
    keyring set https://upload.pypi.org/legacy/ drmjc
* test installation on a fresh osx box

  python -m venv .
  source bin/activate
  bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==0.0.1b16
  
* debugging
  
  bump the version number in _version.py
  ./build.sh
  # wait for test.pypi to index the new package
  # in your venv grab the new version. it'll uninstall the previous one
  source bin/activate
  bin/pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple mity==0.0.1b16
  
# test data

    dx make_download_url --duration 1y kccg-freebayes-mity-resources:/test/kccg-mity-call/inputs/A1.dedup.realigned.recalibrated.chrMT.bam
    dx make_download_url --duration 1y kccg-freebayes-mity-resources:/test/kccg-mity-call/inputs/A1.dedup.realigned.recalibrated.chrMT.bam.bai
    dx make_download_url --duration 1y kccg-freebayes-mity-resources:/assets/hs37d5.fasta-index.tar.gz

# make reference genome
* https://documentation.dnanexus.com/science/scientific-guides/human-genome

    wget https://dl.dnanex.us/F/D/pVG7PjZy4qKBB6ZKbkkF0X6kB0kxf7ZzjpK7fXjY/hs37d5.fasta-index.tar.gz
    tar -xzvf hs37d5.fasta-index.tar.gz; mv genome.dict hs37d5.dict; mv genome.fa hs37d5.fa; mv genome.fa.fai hs37d5.fa.fai
    samtools faidx hs37d5.fa MT -o hs37d5.MT.fa
    samtools faidx hs37d5.MT.fa
    dx upload hs37d5.MT.fa* --path kccg-freebayes-mity-resources:/assets/

    # hg38 from Broad: https://software.broadinstitute.org/gatk/download/bundle
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.dict
    wget https://storage.googleapis.com/genomics-public-data/resources/broad/hg38/v0/Homo_sapiens_assembly38.fasta
    samtools faidx Homo_sapiens_assembly38.fasta 

    samtools faidx Homo_sapiens_assembly38.fasta chrM -o hg38.chrM.fa
    samtools faidx hg38.chrM.fa 

# debug via DNAnexus
```
dx run cloud_workstation --ssh \
  -imax_session_length=1h \
  -ifids=file-BzF4G0j0QpzJP9931xpgKvv2 -ifids=file-BzF4G4Q0Qpz73J1BzVpJ41k4 --yes
```

* see INSTALL.md

# Get testing data

Using the Ashkenazim Trio data.

Getting BAMs:

Links from ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/

```bash
cd test_in
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam.bai

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.hs37d5.2x250.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG003_NA24149_father/NIST_Illumina_2x250bps/novoalign_bams/HG003.hs37d5.2x250.bam.bai

wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.hs37d5.2x250.bam
wget ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG004_NA24143_mother/NIST_Illumina_2x250bps/novoalign_bams/HG004.hs37d5.2x250.bam.bai

samtools view -bh -o HG002.hs37d5.2x250.small.MT.bam HG002.hs37d5.2x250.bam MT:1-700
samtools view -bh -o HG003.hs37d5.2x250.small.MT.bam HG003.hs37d5.2x250.bam MT:1-700
samtools view -bh -o HG004.hs37d5.2x250.small.MT.bam HG004.hs37d5.2x250.bam MT:1-700

samtools index HG002.hs37d5.2x250.small.MT.bam
samtools index HG003.hs37d5.2x250.small.MT.bam
samtools index HG004.hs37d5.2x250.small.MT.bam

```
And need to add read group (RG) for freebayes to get the sample name. See https://www.biostars.org/p/349213/

```bash
cd test_in
samtools addreplacerg -r ID:HG002 -r SM:HG002 HG002.hs37d5.2x250.small.MT.bam -o HG002.hs37d5.2x250.small.MT.RG.bam
samtools addreplacerg -r ID:HG003 -r SM:HG003 HG003.hs37d5.2x250.small.MT.bam -o HG003.hs37d5.2x250.small.MT.RG.bam
samtools addreplacerg -r ID:HG004 -r SM:HG004 HG004.hs37d5.2x250.small.MT.bam -o HG004.hs37d5.2x250.small.MT.RG.bam

samtools index HG002.hs37d5.2x250.small.MT.RG.bam
samtools index HG003.hs37d5.2x250.small.MT.RG.bam
samtools index HG004.hs37d5.2x250.small.MT.RG.bam
```
# Testing the combine variants function in mity normalise
To test the entire function, you need to run mity in a sensitive mode so that there is at least one variant that is repeated.
```bash
mity call --prefix ashkenazim --out-folder-path test_out --min-alternate-fraction 0.001 --region MT:300-500 --normalise test_in/HG002.hs37d5.2x250.small.MT.RG.bam test_in/HG003.hs37d5.2x250.small.MT.RG.bam test_in/HG004.hs37d5.2x250.small.MT.RG.bam 
```

mity call --prefix ashkenazim --out-folder-path test_out --min-alternate-fraction 0.001 --region MT:300-500 --normalise test_in/HG002.hs37d5.2x250.MT.bam

Gives the warning
/Users/putticc/Projects/mity/mitylib/normalise.py:563: RuntimeWarning: divide by zero encountered in log10
  q = round(abs(-10 * log10(1 - binom.cdf(float(AO), float(DP), p))), 2)


# Test GRCh38 from 1000 genomes

wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.cram.crai
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000_genomes_project/data/CEU/NA12878/alignment/NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.bam.bas
wget ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa
samtools view -b -o NA12878.alt_bwamem_GRCh38DH.20150718.CEU.low_coverage.chrM.bam *cram chrM:1-16569

# Docker
```
version=0.4.0
docker build --tag=latest --tag=$version --tag=drmjc/mity:latest --tag=drmjc/mity:$version .
docker push drmjc/mity            # equivalent to docker push drmjc/mity:latest
docker push drmjc/mity:$version

docker run drmjc/mity version
docker run drmjc/mity:latest version
docker run drmjc/mity:$version version

# test mity via Docker
docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity call \
  --prefix ashkenazim \
  --out-folder-path test_out2 \
  --region MT:1-500 \
  --normalise \
  --reference hs37d5 \
  test_in/HG002.hs37d5.2x250.small.MT.RG.bam \
  test_in/HG003.hs37d5.2x250.small.MT.RG.bam \
  test_in/HG004.hs37d5.2x250.small.MT.RG.bam 

docker run -w "$PWD" -v "$PWD":"$PWD" drmjc/mity report \
  --prefix ashkenazim \
  --min_vaf 0.01 \
  --out-folder-path test_out2 \
  test_out2/ashkenazim.mity.vcf.gz
```

# Triple check that joint-calling agrees with single calling
* all agree

    mity call --prefix ashkenazim002 --out-folder-path test_out --min-alternate-fraction 0.001 --region MT:300-320 --normalise test_in/HG002.hs37d5.2x250.small.MT.RG.bam 
    mity call --prefix ashkenazim003 --out-folder-path test_out --min-alternate-fraction 0.001 --region MT:300-320 --normalise test_in/HG003.hs37d5.2x250.small.MT.RG.bam 
    mity call --prefix ashkenazim004 --out-folder-path test_out --min-alternate-fraction 0.001 --region MT:300-320 --normalise test_in/HG004.hs37d5.2x250.small.MT.RG.bam 


    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004   HG003   HG002
    MT      301     .       A       C       1.51592e-12     SBA_FIL DP=51100;MQM=70;MQMR=70;QA=3595;QR=1718745;SAF=18;SAR=137;SRF=24151;SRR=26794;SBR=0.474;SBA=0.116       GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/1:13343:13299,44:13299:442379:33.264:44:1035:23.523:0.0033:31.28      0/1:20895:20851,44:20851:717140:34.394:44:1043:23.705:0.0021:4.82       0/1:16862:16795,67:16795:559226:33.297:67:1517:22.642:0.004:68.7
    MT      302     .       ACCCCCCCTCCCCCGCTTCTGGCCA       CCCCCCCCCTCCCCCCGCTTCTGGCCA     5.09327e+06     POS_FIL DP=107;MQM=70;MQMR=70;QA=2148;QR=509;SAF=0;SAR=93;SRF=3;SRR=11;TYPE=complex;SBR=0.214;SBA=0.0   GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      1/1:47:0,47:0:0:0:47:1064:22.638:1.0:10000      0/0:15:14,1:14:509:36.357:1:21:21.0:0.0667:33.84        1/1:45:0,45:0:0:0:45:1063:23.622:1.0:10000
    MT      302     .       ACCCCCCCTCCCCCGCTTCTGGCCA       CCCCCCCCTCCCCCCGCTTCTGGCCA      5.09327e+06     POS_FIL DP=31;MQM=70;MQMR=70;QA=430;QR=509;SAF=0;SAR=17;SRF=3;SRR=11;TYPE=complex;SBR=0.214;SBA=0.0     GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/0:1:0,1:0:0:0:1:27:27.0:1.0:10000     0/1:27:14,13:14:509:36.357:13:319:24.538:0.4815:159.55  0/0:3:0,3:0:0:0:3:84:28.0:1.0:10000
    
    
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG002
    MT      301     .       A       C       6.07527e-13     SBA_FIL DP=16862;MQM=70;MQMR=70;QA=1517;QR=559226;SAF=9;SAR=58;SRF=7094;SRR=9701;SBR=0.422;SBA=0.134    GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/1:16862:16795,67:16795:559226:33.297:67:1517:22.642:0.004:68.7
    MT      302     .       ACCCCCCCTCCCCCGCTTCTG   CCCCCCCCCTCCCCCCGCTTCTG 1.23996e+06     SBA_FIL;POS_FIL DP=45;MQM=70;MQMR=0;QA=1063;QR=0;SAF=0;SAR=45;SRF=0;SRR=0;TYPE=complex;SBR=0;SBA=0.0    GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      1/1:45:0,45:0:0:0:45:1063:23.622:1.0:10000
    
    
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG003
    MT      301     .       A       C       0       SBA_FIL DP=20895;MQM=70;MQMR=70;QA=1043;QR=717140;SAF=5;SAR=39;SRF=11228;SRR=9623;SBR=0.538;SBA=0.114   GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/1:20895:20851,44:20851:717140:34.394:44:1043:23.705:0.0021:4.82
    MT      302     .       ACCCCCCCTCCCCCGCTTCTGGCCA       CCCCCCCCTCCCCCCGCTTCTGGCCA      525055  POS_FIL DP=27;MQM=70;MQMR=70;QA=319;QR=509;SAF=0;SAR=13;SRF=3;SRR=11;TYPE=complex;SBR=0.214;SBA=0.0     GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/1:27:14,13:14:509:36.357:13:319:24.538:0.4815:159.55
    
    
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  HG004
    MT      301     .       A       C       0       SBA_FIL DP=13343;MQM=70;MQMR=70;QA=1035;QR=442379;SAF=4;SAR=40;SRF=5829;SRR=7470;SBR=0.438;SBA=0.091    GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      0/1:13343:13299,44:13299:442379:33.264:44:1035:23.523:0.0033:31.28
    MT      302     .       ACCCCCCCTCCCCCGCTTCTG   CCCCCCCCCTCCCCCCGCTTCTG 1.10729e+06     SBA_FIL;POS_FIL DP=47;MQM=70;MQMR=0;QA=1064;QR=0;SAF=0;SAR=47;SRF=0;SRR=0;TYPE=complex;SBR=0;SBA=0.0    GT:DP:AD:RO:QR:AQR:AO:QA:AQA:VAF:q      1/1:47:0,47:0:0:0:47:1064:22.638:1.0:10000

# setup dev environment

    pip install --upgrade pip
    # install pip-compile
    pip install wheel
    pip install pip-tools
    pip install twine
    brew install twine
    pip install -r requirements.txt

    python3 -m venv env
    source env/bin/activate
    pip install -r requirements.txt

    # store your twine password in the keychain
    keyring set https://test.pypi.org/legacy/ drmjc
    keyring set https://upload.pypi.org/legacy/ drmjc

# update version in a few places
* mitylib/_version.py
* Docker section above
* Dockerfile
