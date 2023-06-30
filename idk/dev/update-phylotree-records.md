# Phylotree variants
Some of the rows in the mity report were duplicated, and they seemed to happen when there was >1 haplotype record.
Thus I inspected the chr-pos-ref-alt records in the haplotype_data.csv & some were duplicated.
These mostly related to ins records: "Insertions are indicated by the position number preceding the insertion followed by a dot (.), the relative insert position, and the inserted base(s), e.g. "2156.1A".	[phylotree.org]"
Some records were ambiguous.
Looking at all of the rows containing a dot (eg .1T, .2T, .1C etc; n=47) many could be easily fixed. eg 42.1G = 42 T->TG
Some needed manual work: eg 960.XC, 8289.18bpINS (see section headings below)
All the of deletion variants looked ok: "Deletions are indicated by the letter "d" following the (range of) position number(s) involved, e.g. "A249d" or "8281-8289d"."
**Goal**: update the ALT, and where necessary, the POS for these records.

**Methods**:
1. check that the POS-REF was right for the simple to convert variants (eg .1T etc): `samtools mpileup test_in/HG002.hs37d5.2x250.small.MT.RG.bam -r MT:1-17000 -f mitylib/reference/hs37d5.MT.fa | cut -f 1-3 > mitobase.tsv`
2. for each variant:
2.1 google "pyholtree %s" where %s is the phylotree name for the record
2.2 extract the reference sequence. eg `faidx mitylib/reference/hs37d5.MT.fa MT:8261-8301`
2.3 open two of the genbank records hyperlinked in the phylotree page
2.3.1 extract the sequence surrounding the mutation for each record
2.3.2 paste here
2.4 perform manual sequence alignment of the ref sequence to find the location of the insertion
2.5 transform this to hgvs (right-aligned) and vcf-style coordinates (left-aligned)
3. check the deletion variants. the POS needed to be udpated on those with a range.
4. check dupliacates
5. this made me wonder whether the variants needed to be left-aligned as we will be matching vcf records

# Results
## 960.XC
* https://www.phylotree.org/tree/U.htm
>MT:951-981
     GATCACCCC C*TCCCCAAT AAAGCTAAAA CTC
951 agatcacccc cctccccaat aaagctaaaa ctcacctgag ttgtaaaaaa ctccagttga >https://www.ncbi.nlm.nih.gov/nuccore/GU296624
951 agatcacccc cctccccaat aaagctaaaa ctcacctgag ttgtaaaaaa ctccagttga >https://www.ncbi.nlm.nih.gov/nuccore/JQ702746
= 961m.insC
= 960 C CC
= This is thus a clash with phylotree's "960.1C"
>MT:951-981
     GATCACCCC C*TCCCCAAT AAAGCTAAAA CTC
> https://www.ncbi.nlm.nih.gov/nuccore/KC152555
951 agatcacccc cctccccaat aaagctaaaa ctcacctgag ttgtaaaaaa ctccagttga
* Thus merge these 2 records

## 965.XC
> MT:951-981
      GATC ACCCCC***T CCCCAATAAA GCTAAAACTC
gttttagatc accccccccc ccccaataaa gctaaaactc acctgagttg >https://www.ncbi.nlm.nih.gov/nuccore/EU545417          
   ttagatc accccccccc ccccaataaa gctaaaactc acctgagttg >https://www.ncbi.nlm.nih.gov/nuccore/EU545451
= 961m.T>CCCC

## 573.XC
>MT:561-591
      AAAGACA* CCCCCCACAG TTTATGTAGC TTAC
561 ccaaagacac ccccccacag tttatgtagc ttacctcctc aaagcaatac actgaaaatg >https://www.ncbi.nlm.nih.gov/nuccore/JQ705105
561 ccaaagacac ccccccacag tttatgtagc ttacctcctc aaagcaatac actgaaaatg >https://www.ncbi.nlm.nih.gov/nuccore/AY339536
= 574m.insC
= 573 C CC (right aligned)
= 568 C CC (left aligned)

## 8278.XC
* https://www.phylotree.org/tree/R.htm
>MT:8261-8301
ACCCTATAGCACCCCCTCTACCCCCTCTAGAGCCCACTGTA
8271 cacccccccc cctaccccct ctagagccca ctgtaaagct aacttagcat taacctttta >https://www.ncbi.nlm.nih.gov/nuccore/AP010998
**** no signficant BLAST similarity.
8271 tacacaacac taaaggacga acctgatctc ttatactagt atccttaatc atttttattg >https://www.ncbi.nlm.nih.gov/nuccore/DQ112788
**** no signficant BLAST similarity.
^^^ same for other members of this clade: JN857012.1, JN857012, AY255163
* give up on this one.

## 8289.18bpINS
* easy lookup on the phylotree website: 8289.1CCCCCTCTACCCCCTCTA

## 8289.9bpINS
* easy lookup on the phylotree website: 8289.1CCCCCTCTA

## fix the POS for the deletion records

## check duplicates CHR-POS-ALT
* 13 duplicates...
* all turn out to be the !! records, which are: "Mutations that are reversions to an ancestral state (back mutations) are indicated with an exclamation mark (!), two exclamation marks for a double back mutation (!!), etc., e.g. "A15301G!"."
* in all cases, resolve the duplicate by taking the `phylotree_haplotype` from the !! allele and adding to the `phylotree_haplotype` record for the non !! allele.
$ R
library(mjcbase)
df <- pbpaste("data.frame")
df[which(duplicated(df)),]
     REF ALT   POS  <MJC: The allele listed below will be merged into the non '!!' allele>
55     G   A   103  * G103A!!
76     C   T   146  * C146T!!
83     T   C   152  * T152C!!!
85     C   T   152  * C152T!! (wait how can REF = T and C. check re the reversion to ancestral allele.)
94     C   T   182  * C182T!!
116    C   T   195  * C195T!!
2659   G   A 10398  * G10398A!!
3930   G   A 15301  * G15301A!!
4247   A   G 16129  * A16129G!!
4315   C   T 16189  * C16189T!!
4429   T   C 16278  * T16278C!!
4477   C   T 16311  * C16311T!!

## Turn this into a VCF
    ##fileformat=VCFv4.2							
    ##fileDate=20191106							
    ##reference=/Users/mcowley/src/mity/mitylib/reference/hs37d5.MT.fa							
    ##contig=<ID=MT,length=16569>							
    ##ALT=<ID=NON_REF,Description="Represents any possible alternative allele at this location">							
    ##INFO=<ID=PHYLOTREE,Number=1,Type=String,Description="Phylotree record ID">							
    #CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
    MT	10	T10C	T	C	.	PASS	PHYLOTREE=J1b1a1c
    MT	16	A16t	A	T	.	PASS	PHYLOTREE=G3a1,K1a11
    MT	26	C26T	C	T	.	PASS	PHYLOTREE=R0a1b
    MT	41	T41C!	T	C	.	PASS	PHYLOTREE=T2e1a1a
    MT	41	C41T	C	T	.	PASS	PHYLOTREE=H6a1a2a,T2b34,T2e1
    MT	42	42.1G	T	TG	.	PASS	PHYLOTREE=H1a3a2,N9a10b
    ...
## left-normalise this file
* mity's vcf will have left-normalised variants. Likely that there are at least some records in the Phylotree file that will need to be left-normalised.

    vt normalize -n -r mitylib/reference/hs37d5.MT.fa -o dev/phylotree.cleaned.left.vcf dev/phylotree.cleaned.vcf 

* lots of REF inconsistencies! All appear to be del records, so probably out-by-one errors.
reference bases not consistent: MT:70-71  GG(REF) vs GT(FASTA)                     # 
reference bases not consistent: MT:104-110  CCGGAGC(REF) vs CGGAGCA(FASTA)         # 
reference bases not consistent: MT:105-111  CGGAGCA(REF) vs GGAGCAC(FASTA)         # 
reference bases not consistent: MT:248-249  AA(REF) vs AT(FASTA)                   # 
reference bases not consistent: MT:289-291  AAA(REF) vs AAT(FASTA)                 # 
reference bases not consistent: MT:290-294  AATTT(REF) vs ATTTC(FASTA)             # 
reference bases not consistent: MT:298-299  CC(REF) vs CA(FASTA)                   # 
reference bases not consistent: MT:308-309  CC(REF) vs CT(FASTA)                   # 
reference bases not consistent: MT:336-337  AA(REF) vs AC(FASTA)                   # 
reference bases not consistent: MT:454-455  TT(REF) vs TC(FASTA)                   # 
reference bases not consistent: MT:458-459  CC(REF) vs CT(FASTA)                   # 
reference bases not consistent: MT:497-498  CC(REF) vs CG(FASTA)                   # 
reference bases not consistent: MT:959-960  CC(REF) vs CT(FASTA)                   # 
reference bases not consistent: MT:1408-1409  AA(REF) vs AG(FASTA)                 # 
reference bases not consistent: MT:1655-1656  AA(REF) vs AC(FASTA)                 # 
reference bases not consistent: MT:2073-2074  AA(REF) vs AT(FASTA)                 # 
reference bases not consistent: MT:2394-2395  AA(REF) vs AC(FASTA)                 # 
reference bases not consistent: MT:4316-4317  AA(REF) vs AC(FASTA)                 # 
reference bases not consistent: MT:5751-5752  AA(REF) vs AG(FASTA)                 # 
reference bases not consistent: MT:5893-5894  CA(REF) vs AC(FASTA)                 # 
reference bases not consistent: MT:5898-5899  CC(REF) vs CA(FASTA)                 # 
reference bases not consistent: MT:5959-5961  TTT(REF) vs CCT(FASTA)               # 
reference bases not consistent: MT:6464-6466  TCT(REF) vs GTC(FASTA)               # 
reference bases not consistent: MT:7470-7471  CC(REF) vs CA(FASTA)                 # 
reference bases not consistent: MT:7502-7503  CC(REF) vs CA(FASTA)                 # 
reference bases not consistent: MT:8280-8289  ACCCCCTCTA(REF) vs CCCCCTCTAG(FASTA) # 
reference bases not consistent: MT:15943-15944  TT(REF) vs TC(FASTA)               # 
reference bases not consistent: MT:16165-16166  AA(REF) vs AC(FASTA)               # 
reference bases not consistent: MT:16192-16193  CC(REF) vs CA(FASTA)               # 
reference bases not consistent: MT:16256-16257  CC(REF) vs CA(FASTA)               # 
reference bases not consistent: MT:16324-16325  TT(REF) vs TA(FASTA)               # 

** fix all these

    vt normalize -n -r mitylib/reference/hs37d5.MT.fa -o dev/phylotree.cleaned.left.vcf dev/phylotree.cleaned.vcf 
...
          no. left aligned                      : 65
       total no. variants normalized            : 65
       total no. variants observed              : 4546

** convert back to Clare's phylotree format
grep -v ^# dev/phylotree.cleaned.vcf | awk -v OFS='\t' 'BEGIN{print "phylotree_mut\tREF\tALT\tPOS\tphylotree_haplotype\tCHR"}; {print $3, $4, $5, $2, $8, $1}' | perl -p -e 's/PHYLOTREE=//' > dev/haplotype_data.tsv

python
import pandas as pd 
tsv_file='dev/haplotype_data.tsv'
csv_table=pd.read_table(tsv_file,sep='\t')
csv_table.to_csv('dev/haplotype_data.csv',index=False)
less dev/haplotype_data.csv

** compare format
* almost identical. imporantly the records with commas in them are enclosed with ""

    $ head dev/haplotype_data.csv
    phylotree_mut,REF,ALT,POS,phylotree_haplotype,CHR
    T10C,T,C,10,J1b1a1c,MT
    A16t,A,T,16,"G3a1,K1a11",MT
    C26T,C,T,26,R0a1b,MT
    T41C!,T,C,41,T2e1a1a,MT
    C41T,C,T,41,"H6a1a2a,T2b34,T2e1",MT
    42.1G,T,TG,42,"H1a3a2,N9a10b",MT
    44.1C,C,CC,44,"D5a2a1a,H15a1b,H6b1,L1c1a1a,X2j",MT
    G47A,G,A,47,C4a2c,MT
    G53A,G,A,53,F1a3b,MT

    $ head mitylib/annot/haplotype_data.csv
    phylotree_mut,REF,ALT,POS,phylotree_haplotype,CHR
    A235G,A,G,235,"A,B4a1a1a8,F2d,J2a2a,K3,L2a2b2,L3k,Q1c2a,U6b3a,X2m1",MT
    A663G,A,G,663,"A,U1a1c1c1",MT
    A1736G,A,G,1736,A,MT
    T4248C,T,C,4248,"A,D4a4,E1a",MT
    A4824G,A,G,4824,"A,B5a2a2b1,L1c1b1",MT
    C8794T,C,T,8794,A,MT
    C16290T,C,T,16290,"A,D1g4,D4o,G1a2,H1k,J1b1b1a,J1b9,L0d3,L2a1b1a,L2a1c3a1,M12,M22a,M35c,N13,R23,V1b",MT
    G16319A,G,A,16319,"A,B2c2a,B5b1a2,C1c3,D4b1,D4d,D4o1,G2a4,H13a2c,I1a1a3,I1c,J1c8a,J2b1c1,K1b1a,L0a4,L3d1a1'2,M12b1a2,M2a'b,M32,M35a1a,M3c1b1a,M40a1,M47,M55,M60a,M6a2,M79,M7c1a3,M7c1c1,M8a,N2a,N5a,S4,N1a3a1,P4a,R7,U5b1f1a,X2a1a1",MT
    G1442A,G,A,1442,"A1,A2aj,B4f,K1c1d,L2b'c,L3e1g,M7a1a2,R7a'b",MT

** overwrite and test

    cp mitylib/annot/haplotype_data.csv dev/haplotype_data.original.csv
    cp dev/haplotype_data.csv mitylib/annot/