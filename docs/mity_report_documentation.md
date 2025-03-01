# mity report documentation

## Description of column headings

| Term                                              | Definition and Information                                                                                                                |
| ------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------- |
| SAMPLE                                            | Sample name from the VCF (Variant Call Format) file.                                                                                      |
| HGVS                                              | Variant annotation using HGVS (Human Genome Variation Society) syntax for SNPs, insertions, or deletions.                                 |
| GENE/LOCUS                                        | The gene the variant is located in, or the corresponding MITOMAP entry if not in a gene.                                                  |
| GENE/LOCUS DESCRIPTION                            | Description of the gene or locus where the variant is found.                                                                              |
| TOTAL LOCUS DEPTH                                 | The total sequencing depth (number of reads) at the locus of interest.                                                                    |
| VARIANT HETEROPLASMY                              | The proportion of alternate allele reads to total reads: alt_depth / (ref_depth + alt_depth).                                             |
| ALT DEPTH                                         | The number of sequencing reads that contain the alternate (variant) allele.                                                               |
| REF DEPTH                                         | The number of sequencing reads that contain the reference allele.                                                                         |
| TOTAL SAMPLE DEPTH                                | The total number of sequencing reads for the sample at the given position: alt_depth + ref_depth.                                         |
| VARIANT QUALITY                                   | The Phred-scaled quality score of the variant call, representing confidence in the call.                                                  |
| TIER                                              | Categorization based on heteroplasmy fraction: Tier 1 (VAF ≥ 1%), Tier 2 (VAF < 1% & alt_depth ≥ 10), Tier 3 (VAF < 1% & alt_depth < 10). |
| COHORT COUNT                                      | The number of times the variant (defined by position and alternate allele) appears in all input VCFs.                                     |
| COHORT FREQUENCY                                  | The frequency of the variant in the cohort: COHORT COUNT / TOTAL SAMPLES.                                                                 |
| MQM_INFO                                          | Mean mapping quality of observed alternate alleles.                                                                                       |
| MQMR_INFO                                         | Mean mapping quality for reference alleles.                                                                                               |
| QA_INFO                                           | Sum of Phred-scaled quality scores for alternate allele reads.                                                                            |
| QR_INFO                                           | Sum of Phred-scaled quality scores for reference allele reads.                                                                            |
| SAF_INFO                                          | Number of alternate allele reads observed on the forward strand.                                                                          |
| SAR_INFO                                          | Number of alternate allele reads observed on the reverse strand.                                                                          |
| SRF_INFO                                          | Number of reference allele reads observed on the forward strand.                                                                          |
| SRR_INFO                                          | Number of reference allele reads observed on the reverse strand.                                                                          |
| SBR_INFO                                          | Signal Bias Ratio: For all alleles, if RO > 15 and (SBR > 0.9 or SBR < 0.1), it indicates strand bias.                                    |
| SBA_INFO                                          | Signal Bias for Alternate alleles: If AO > 15 and (SBA > 0.9 or SBA < 0.1), it indicates strand bias.                                     |
| POS_FILTER                                        | Indicates whether the variant falls within a known problematic region, such as MT:302-319 or MT:3105-3108.                                |
| SBR_FILTER                                        | A filter flag for strand bias in reference alleles.                                                                                       |
| SBA_FILTER                                        | A filter flag for strand bias in alternate alleles.                                                                                       |
| MQMR_FILTER                                       | A filter flag applied when MQMR (reference mapping quality) is below 30.                                                                  |
| AQR_FILTER                                        | A filter flag applied when AQR (average quality of reference reads) is below 20.                                                          |
| ANTICODON                                         | The anticodon sequence if the variant is located in a tRNA gene.                                                                          |
| GENE                                              | The gene name where the variant occurs.                                                                                                   |
| GENE BIOTYPE                                      | The biotype classification of the gene (e.g., protein-coding, tRNA, rRNA).                                                                |
| PHYLOTREE MUT                                     | The mutation annotation from the PhyloTree database.                                                                                      |
| PHYLOTREE HAPLOTYPE                               | The associated haplogroup assignment based on PhyloTree.                                                                                  |
| MGRB FILTER                                       | Indicates whether the variant passes filtering criteria in the MGRB dataset.                                                              |
| MGRB AN                                           | Allele number in the MGRB dataset.                                                                                                        |
| MGRB AC                                           | Allele count in the MGRB dataset.                                                                                                         |
| MGRB FREQUENCY                                    | Frequency of the variant allele in the MGRB dataset.                                                                                      |
| MITOMAP DISEASE AC                                | Allele count of the variant in MITOMAP’s disease-associated dataset.                                                                      |
| MITOMAP DISEASE AF                                | Allele frequency of the variant in MITOMAP’s disease-associated dataset.                                                                  |
| MITOMAP DISEASE AACHANGE                          | The associated amino acid change if applicable.                                                                                           |
| MITOMAP DISEASE HOMOPLASMY                        | Whether the variant has been observed in a homoplasmic state in disease cases.                                                            |
| MITOMAP DISEASE HETEROPLASMY                      | Whether the variant has been observed in a heteroplasmic state in disease cases.                                                          |
| MITOMAP DISEASE PUBMED IDS                        | List of PubMed IDs of studies linking the variant to disease.                                                                             |
| MITOMAP DISEASE DISEASE                           | Name of the disease(s) associated with the variant.                                                                                       |
| MITOMAP DISEASE DISEASE STATUS                    | Classification of the disease association (e.g., confirmed, probable, uncertain).                                                         |
| MITOMAP DISEASE HGFL                              | The Human Gene and Disease Feature List (HGFL) entry for mitochondrial diseases associated with specific variants.                        |
| MITOMAP CONFIRMED MUTATIONS LOCUS                 | The specific mitochondrial gene or locus where a confirmed pathogenic mutation occurs.                                                    |
| MITOMAP CONFIRMED MUTATIONS LOCUSTYPE             | The classification of the locus (e.g., protein-coding gene, rRNA, tRNA, or control region).                                               |
| MITOMAP CONFIRMED MUTATIONS ASSOCIATEDDISEASE     | The disease or condition associated with the confirmed mitochondrial mutation.                                                            |
| MITOMAP CONFIRMED MUTATIONS ALLELE                | The specific allele (variant) that has been confirmed as pathogenic in mitochondrial disease.                                             |
| MITOMAP CONFIRMED MUTATIONS AMINOACIDCHANGE       | The change in the amino acid sequence resulting from a confirmed pathogenic mutation, if applicable.                                      |
| MITOMAP CONFIRMED MUTATIONS STATUSMITOMAPCLINGEN  | The classification status of the mutation in MITOMAP and ClinGen, indicating whether it is pathogenic, likely pathogenic, or uncertain.   |
| MITOMAP CONFIRMED MUTATIONS LASTUPDATE            | The date when the mutation entry was last updated in MITOMAP.                                                                             |
| MITOMAP MUTATIONS CODING CONTROL LOCUS            | The mitochondrial gene or control region affected by a mutation.                                                                          |
| MITOMAP MUTATIONS CODING CONTROL ALLELE           | The specific allele involved in the mutation within the coding or control region.                                                         |
| MITOMAP MUTATIONS CODING CONTROL DISEASE          | The disease or disorder linked to the mutation occurring in the coding or control region.                                                 |
| MITOMAP MUTATIONS CODING CONTROL NUCLEOTIDECHANGE | The specific nucleotide-level alteration in the coding or control region of the mitochondrial genome.                                     |
| MITOMAP MUTATIONS CODING CONTROL AMINOACIDCHANGE  | The resulting change in the amino acid sequence due to a coding-region mutation, if applicable.                                           |
| MITOMAP MUTATIONS CODING CONTROL PLASMY           | The type of heteroplasmy (heteroplasmic or homoplasmic) observed for the mutation.                                                        |
| MITOMAP MUTATIONS CODING CONTROL STATUS           | The classification status of the mutation (e.g., pathogenic, likely pathogenic, uncertain, or benign).                                    |
| MITOMAP MUTATIONS CODING CONTROL GB FREQ          | The frequency of the mutation in GenBank mitochondrial sequences.                                                                         |
| MITOMAP MUTATIONS CODING CONTROL GB SEQS          | The number of GenBank mitochondrial genome sequences in which the mutation has been observed.                                             |
| MITOMAP MUTATIONS CODING CONTROL REFERENCES       | Research articles or curated references that discuss the mutation's significance.                                                         |
| MITOMAP MUTATIONS RNA LOCUS                       | The specific mitochondrial RNA gene (tRNA or rRNA) affected by a mutation.                                                                |
| MITOMAP MUTATIONS RNA DISEASE                     | The disease or disorder associated with the mitochondrial RNA mutation.                                                                   |
| MITOMAP MUTATIONS RNA ALLELE                      | The specific allele variant of the mitochondrial RNA gene that has been mutated.                                                          |
| MITOMAP MUTATIONS RNA RNA                         | The affected RNA sequence or structural change resulting from the mutation.                                                               |
| MITOMAP MUTATIONS RNA HOMOPLASMY                  | Indicates whether the mutation is present in all copies of mitochondrial DNA in a cell (homoplasmic state).                               |
| MITOMAP MUTATIONS RNA HETEROPLASMY                | Indicates whether the mutation is present in only a subset of mitochondrial DNA copies in a cell (heteroplasmic state).                   |
| MITOMAP MUTATIONS RNA STATUS                      | The classification status of the RNA mutation (e.g., pathogenic, likely pathogenic, uncertain).                                           |
| MITOMAP MUTATIONS RNA MITOTIP                     | MitoTIP score, a computational prediction tool used to assess the pathogenicity of mitochondrial tRNA mutations.                          |
| MITOMAP MUTATIONS RNA GB FREQ                     | The frequency of the RNA mutation in GenBank mitochondrial sequences.                                                                     |
| MITOMAP MUTATIONS RNA GB SEQS                     | The number of GenBank mitochondrial genome sequences in which the RNA mutation has been observed.                                         |
| MITOMAP MUTATIONS RNA REFERENCES                  | Published studies or curated references discussing the mitochondrial RNA mutation and its biological significance.                        |
| MITOMAP POLYMORPHISMS AC                          | Allele count in MITOMAP’s polymorphism dataset.                                                                                           |
| MITOMAP POLYMORPHISMS AF                          | Allele frequency in MITOMAP’s polymorphism dataset.                                                                                       |
| MITOMAP POLYMORPHISMS HGFL                        | Haplogroup Frequency List, indicating common haplogroups for the variant.                                                                 |
| MITOMAP VARIANTS CODING LOCUS                     | The specific gene or locus in the mitochondrial genome where a coding-region variant occurs.                                              |
| MITOMAP VARIANTS CODING NUCLEOTIDECHANGE          | The specific nucleotide alteration (substitution, deletion, or insertion) observed in the coding region of the mitochondrial genome.      |
| MITOMAP VARIANTS CODING CODONNUMBER               | The codon number in the gene where the variant occurs, indicating its position in the coding sequence.                                    |
| MITOMAP VARIANTS CODING CODONPOSITION             | The specific position within the codon affected by the nucleotide change (first, second, or third base of the codon).                     |
| MITOMAP VARIANTS CODING AMINOACIDCHANGE           | The resulting amino acid substitution due to the nucleotide change, if it leads to a missense or nonsense mutation.                       |
| MITOMAP VARIANTS CODING GB FREQ                   | The frequency of the variant observed in GenBank mitochondrial genome sequences, indicating how common it is in public databases.         |
| MITOMAP VARIANTS CODING GB SEQS                   | The number of GenBank mitochondrial genome sequences in which the variant has been observed.                                              |
| MITOMAP VARIANTS CODING CURATEDREFERENCES         | Published research articles or curated references that discuss the variant in the coding region, supporting its biological significance.  |
| MITOMAP VARIANTS CONTROL LOCUS                    | The specific locus in the mitochondrial control region where a variant occurs.                                                            |
| MITOMAP VARIANTS CONTROL NUCLEOTIDECHANGE         | The nucleotide change observed in the control region of the mitochondrial genome.                                                         |
| MITOMAP VARIANTS CONTROL GB FREQ                  | The frequency of the control-region variant in GenBank mitochondrial genome sequences.                                                    |
| MITOMAP VARIANTS CONTROL GB SEQS                  | The number of GenBank sequences containing this variant in the control region.                                                            |
| MITOMAP VARIANTS CONTROL CURATED REFERENCES       | Research articles or curated sources discussing the significance of the control-region variant.                                           |
| MITOTIP SCORE                                     | A computational pathogenicity score for tRNA variants.                                                                                    |
| MITOTIP PERCENTILE                                | Percentile ranking of the MITOTIP score.                                                                                                  |
| MITOTIP QUARTILE                                  | Quartile-based ranking of MITOTIP score (Q1-Q4).                                                                                          |
| MITOTIP SCORE INTERPRETATION                      | Interpretation of the pathogenicity score based on its percentile.                                                                        |
| MITOMAP STATUS                                    | The strength of evidence supporting the MITOMAP annotation.                                                                               |
| COUNT                                             | Count of MITOTIP-associated observations.                                                                                                 |
| PERCENTAGE                                        | Percentage of MITOTIP observations.                                                                                                       |
| GT_FORMAT                                         | Genotype format field in the VCF.                                                                                                         |
| QR_FORMAT                                         | Sum of Phred-scaled base quality scores for reference allele reads (per sample).                                                          |
| AQR_FORMAT                                        | Average base quality of the reference reads: AQR = QR / RO.                                                                               |
| QA_FORMAT                                         | Sum of Phred-scaled base quality scores for alternate allele reads (per sample).                                                          |
| AQA_FORMAT                                        | Average base quality of the alternate reads: AQA = QA / AO.                                                                               |

## Mitomap source links and conversions

| Mitomap Annotation File                 | Source                                                                  | Python script to generate vcf                      |
| --------------------------------------- | ----------------------------------------------------------------------- | -------------------------------------------------- |
| mitomap_disease.vcf.gz                  | https://www.mitomap.org/foswiki/bin/view/MITOMAP/Resources              |                                                    |
| mitomap_polymorphisms.vcf.gz            | https://www.mitomap.org/foswiki/bin/view/MITOMAP/Resources              |                                                    |
| mitotip_scores.vcf.gz                   | https://www.mitomap.org/foswiki/bin/view/MITOMAP/MitoTipInfo            | convert_mitomap_scores_to_vcf.py                   |
| mitomap_mutations_rna.vcf.gz            | https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsRNA           | convert_mitomap_mutations_rna_to_vcf.py            |
| mitomap_mutations_coding_control.vcf.gz | https://www.mitomap.org/foswiki/bin/view/MITOMAP/MutationsCodingControl | convert_mitomap_mutations_coding_control_to_vcf.py |
| mitomap_variants_control.vcf.gz         | https://www.mitomap.org/foswiki/bin/view/MITOMAP/VariantsControl        | convert_mitomap_variants_control_to_vcf.py         |
| mitomap_variants_coding.vcf.gz          | https://www.mitomap.org/foswiki/bin/view/MITOMAP/VariantsCoding         | convert_mitomap_variants_coding_to_vcf.py          |
| mitomap_confirmed_mutations.vcf.gz      | https://www.mitomap.org/foswiki/bin/view/MITOMAP/ConfirmedMutations     | convert_mitomap_confirmed_mutations_to_vcf.py      |

Most of the mitomap sources come in other forms, e.g. `csv` or `tsv`, but we need `vcf`s to properly annotate with `vcfanno` so in `tools/` we have a series of conversion scripts unique to each mitomap source.

### Notes

`mitomap_variants_control.vcf.gz` from the original `VariantsControl MITOMAP Foswiki.csv` has 4 lines which were somewhat difficult to convert to the vcf version. Namely these lines:

```csv
"568","MT-HV3","C-C(2-8)","0.000%(0.000%)","0","2"
"573","MT-HV3","C-C(2-8)","0.000%(0.000%)","0","40"
"16184","MT-HV1","C-C(2-5)","0.000%(0.000%)","0","12"
"16193","MT-HV1","C-C(2-3)","0.000%(0.000%)","0","18"
```
