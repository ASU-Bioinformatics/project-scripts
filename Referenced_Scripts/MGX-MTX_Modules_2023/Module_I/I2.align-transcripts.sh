#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 0-12:00
#SBATCH -c 1
#SBATCH --mem=64G

#!!! this will need to be changed to HISAT2, since STAR didn't work
#!!! and minimap2 says it isn't designed for RNA reads
#!!! which I believe since all the alignments ended up as 0s...

#!!! for strandedness option for HISAT2, I think the KAPA kit has read1 as antisense and read2 as sense, like TruSeq


#### set up environment ####

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

#### read in variables ####

transcriptDir=$1
hisatIndex=$2
sampleID=$3
outputDir=$4

cd $outputDir
/data/biocore/programs/hisat2-2.2.1/hisat2 --dta -p 64 \
                                           --summary-file "$outputDir"/"$sampleID-summary.txt" \
                                           -x "$hisatIndex" \
                                           -1 "$transcriptDir"/"$sampleID"_SQP_L001_R1_001.fastq.gz \
                                           -2 "$transcriptDir"/"$sampleID"_SQP_L001_R2_001.fastq.gz \
                                           -S "$outputDir"/"$sampleID-transcripts-on-genes.sam"

samtools view -bS "$outputDir"/"$sampleID-transcripts-on-genes.sam" > "$outputDir"/"$sampleID-transcripts-on-genes.bam"

samtools sort "$outputDir"/"$sampleID-transcripts-on-genes.bam" -o "$outputDir"/"$sampleID-transcripts-on-genes.sorted.bam"

samtools index "$outputDir"/"$sampleID-transcripts-on-genes.sorted.bam"
