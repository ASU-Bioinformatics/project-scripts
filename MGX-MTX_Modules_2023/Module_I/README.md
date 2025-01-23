# Module I

This module aligns MTX data to the MGX coassembly, creates a gene count matrix for the sample set, and evaluates those transcripts for differential expression and functional enrichment using a standard RNAseq pipeline (DESeq2, edgeR, NOISeq, clusterProfiler, etc.). The first half of the module is primarily in bash, and the second half is primarily in R. Because the R scripts can have high memory demands, they have an accompanying wrapper script that can run the R programs through the Slurm scheduler on Sol.

## I1.hisat-index.sh

This script uses the RNA aligner [HiSat2](https://daehwankimlab.github.io/hisat2/,"Title") to create an index for the MGX coassembly fasta file. To ensure downstream compatibility with the GFF file and HTSeq counts, the contig names for the coassembly are truncated at the first whitespace. A large index is automatically generated as the coassembly fasta file typically contains a significant number of contigs and the default attempt to create a small index will lead to an error.

On Sol, this module relies on the Conda environment `/data/biocore/programs/mamba-envs/biocore-rna`; if you require a local installation, the included file `biocoreRNA.yml` can be used to create your own environment. HiSat2 itself is not installed as a part of the environment and must be separately installed. If you have access to the Biocore drive there is an installation available at `/data/biocore/programs/hisat2-2.2.1/`.

A sample call:

    sbatch ./I1.hisat-index.sh \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-index \
      --prefix final.contigs.reference \
      --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa

Input parameters:
* `-o,--outDir`: The name of the directory to which the shortname contig fasta file and HiSat2 index should be written
* `-p,--prefix`: The prefix for the HiSat2 index (this is how the index will be referenced in downstream scripts as well)
* `-c,--contigs`: The pathname for the MGX coassembly fasta to be indexed

Output files include a large HiSat2 index for the MGX coassembly with the specified prefix, and the copy of the coassembly fasta with truncated contig names.

## I2.align-transcripts.wrapper.sh

This script uses HiSat2 to align the MTX transcript reads to the MGX coassembly, then converts the resulting SAM files to BAM format and sorts the BAM files to prepare them for abundance profiling in the next step. To increase the efficiency of the alignment step, a wrapper file is used to generate sbatch calls for each sample specified in a directory or list.

On Sol, this module relies on the Conda environment `/data/biocore/programs/mamba-envs/biocore-rna`; if you require a local installation, the included file `biocoreRNA.yml` can be used to create your own environment. HiSat2 itself is not installed as a part of the environment and must be separately installed. If you have access to the Biocore drive there is an installation available at `/data/biocore/programs/hisat2-2.2.1/`.

A sample call:

    bash ./I2.align-transcripts.wrapper.sh \
      --transcriptDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-alignments \
      --hisatIndex /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-index/final.contigs.reference

Input parameters:
* `-t,--transcriptDir`: the directory containing fastq.gz files for MTX reads, assuming paired-end data.
* `-o,--outDir`: the directory to which the alignment SAM, BAM, and summary stats files should be written
* `-h,--hisatIndex`: the path name to the HiSat2 index, including the index prefix
* `-l,--list`: a list of sample IDs to run alignment for. Sample IDs should include the full name of the fastq.gz file up to `_SQP_L001_R*_001.fastq.gz` in the standard naming convention.

Output files include a sam, bam, sorted bam, bam.bai index, and summary text file for each sample. The summary file provides statistics for the number of reads aligned paired and singly as well as the overall alignment rate.

## I3.abundance-profiling.wrapper.sh

This script uses [HTSeq count](https://htseq.readthedocs.io/en/master/htseqcount.html,"Title") to determine the number of reads in each sample that are aligned to each gene and output it in tabular format - essentially, an abundance count of the alignment data stored in the sorted BAM files from the previous script.

On Sol this can be run with the Conda environment `/data/biocore/programs/mamba-envs/htseq-env` if you have access to the Biocore data drive; if not, you can create a local installation using the included `htseq.yml` file.

A sample call:

    bash ./I3.abundance-profiling.wrapper.sh \
      --alignmentDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-alignments \
      --gffFile /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/out.CAT.predicted_proteins.gff \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts

Input parameters:
* `-a,--alignmentDir`: the directory containing the sorted BAM files and accompanying BAI index files
* `-g,--gffFile`: the path name to a GFF file containing the predicted protein annotations on the coassembled contigs
* `-o,--outDir`: the directory to which the count files should be written
* `-l,--list`: a list of sample IDs within the alignment directory that should be counted, if only some of the samples in the directory need counts

Output files include a tab-delimited transcript counts file for each sample.

## I4.merge-htseq.sh

This script takes the individual counts files from the previous script and merges them to create a unified gene count matrix for differential expression analysis. It is designed to merge all files in the specified directory that have the suffix `.htseqcounts.txt`, which is the automatic suffix for the previous script. Since only one input argument is needed, it is positional rather than flagged.

On Sol this can be run with the Conda environment `/data/biocore/programs/mamba-envs/htseq-env` if you have access to the Biocore data drive; if not, you can create a local installation using the included `htseq.yml` file.

A sample call:

    sbatch ./I4.merge-htseq.sh /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts

The input is the directory containing the HTseq counts files for the individual samples in the data set, and the output is a tab-delimited count matrix.

## Before moving on
Check the alignment rates in the summary text files for each sample. High alignment rates (>95%) show that the coassembly file was able to incorporate the genetic diversity of the corresponding MGX sample, and that the transcripts present in the MTX sample do not contain reads not represented by the MGX sample. Low alignment rates are troubling because they suggest many transcripts whose corresponding genes were not present in the MGX data, which doesn't make biological sense and could indicate contamination or degradation.
