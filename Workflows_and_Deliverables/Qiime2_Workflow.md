*Last updated September 1, 2026*

## Analysis of 16S/18S/ITS Amplicon Regions from Illumina Short-Read Sequencing

This document describes the standard workflow used in ASU's bioinformatics core facility for basic analysis of short-read amplicon sequencing on standard taxonomically-relevant genomic regions.

### Sequencing Input and QC

This section of the analysis is performed for all amplicon sequencing, even if no bioinformatics is requested, to ensure that the returned fastq reads have the expected read depth and quality and can be used without manually removing adapter sequences.

#### Evaluate Raw Reads

FastQC (v0.12.1) and MultiQC (v1.20) are used to ensure the raw reads have sufficient depth and acceptable quality. Because the adapter removal step will reduce the overall read depth and length with poor quality samples, at this step we are primarily looking for samples where sequencing has clearly failed or has fewer than 10K reads. These samples will be rerun at the Desert Southwest Genome Center (DSGC).

#### Remove Adapter Sequence

Because the method of library prep performed at the DSGC for amplicon regions leaves a 5' adapter sequence at the beginning of each read in addition to the more common read-through adapter artifact at the 3' end of shorter inserts, we perform adapter removal with cutadapt (v4.8) using the following sequences and a minimum retained length of 10bp:

| amplicon region | sequence                                                            |
| --------------- | ------------------------------------------------------------------- |
| 3' all regions  | CTGTCTCTTATA                                                        |
| 5' 16S R1       | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGGTGCCAGCMGCCGCGGTAA                |
| 5' 16S R2       | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGGGACTACHVGGGTWTCTAAT              |
| 5' 18S NOS R1   | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTATCGCCGTTCGGTACACACCGCCCGTC       |
| 5' 18S NOS R2   | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTCAGTCAGCATGATCCTTCTGCAGGTT     |
| 5' 18S v4 R1    | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTATCGCCGTTCGCCAGCASCYGCGGTAATTCC   |
| 5' 18S v4 R2    | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTCAGTCAGCAACTTTCGTTCTTGATYRA    |
| 5' 18S v9 R1    | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTATCGCCGTTCGTTGTACACACCGCCC        |
| 5' 18S v9 R2    | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTCAGTCAGCACCTTCYGCAGGTTCACCTAC  |
| 5' ITS R1       | TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTCTTGGTCATTTAGAGGAAGTAA           |
| 5' ITS R2       | GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGATGCTGCGTTCTTCATCGATGC            |

After adapter removal we rerun FastQC and MultiQC to verify that all samples still have more than 10K reads (or another read depth requested prior to sequencing) without significant dips in quality over the course of the reads.

If specifically requested (for example, to match maximum read length with previous runs from shorter-read kits), we will also cut the reads following adapter removal to the desired length.

### Qiime2

Secondary analysis of conserved region amplicon sequencing is carried out using the Qiime2 (v2025.7) command line interface.

#### Interleave and Denoise Reads

Since the DSGC provides paired-end sequences, we use the `denoise-paired` function within the Qiime2 `dada2` module to interleave the forward and reverse reads, denoise, and remove chimeras. Using the `feature-table` module, we again evaluate the sample quality and read depth. If the read depth is below 10K following the dada2 QC filter, we ask the DSGC to rerun the samples to increase the number of reads continuing into the analysis.

#### Calculate Alpha and Beta Diversity

The Qiime2 module `phylogeny` is used to generate a single rooted tree for the experiment so that weighted beta diversity metrics can be calculated. This tree is built using the individual ASVs from each sample, so it is not dependent on taxonomic classifications.

Beta and alpha diversity are both calculated using the `diversity` module, specifically including the Jaccard, weighted and unweighted UniFrac, and Bray-Curtis metrics for beta diversity as well as Faith's PD, evenness, and the `diversity` module's alpha-rarefaction command (which includes feature count and the Shannon index) for alpha diversity. All the beta diversity metrics are graphed as 3D interactive emperor plots as well as with pairwise comparison between metadata groups. The alpha rarefaction output is graphed as a curve to evaluate whether the read depth is sufficient to characterize the population.
