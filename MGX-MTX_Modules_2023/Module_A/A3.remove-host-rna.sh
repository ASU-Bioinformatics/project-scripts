#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.%a.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.%a.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 2-0:00
#SBATCH -c 6
#SBATCH --mem=64G

### Define Variables ####

sid=$1
inputDir=$2
refDir=$3
tempDir=$4
outputDir=$5

cd "$inputDir"
cp "$sid"_*.fastq.gz "$tempDir"

cd "$tempDir"
STAR \
    --genomeDir "$refDir" \
    --readFilesCommand gunzip -c \
    --readFilesIn "$inputDir"/"$sid"_SQP_L001_R1_001.fastq.gz "$inputDir"/"$sid"_SQP_L001_R2_001.fastq.gz \
    --outFileNamePrefix "$sid"_STAR_ \
    --runThreadN 12 \
    --outSAMtype BAM SortedByCoordinate \
    --outSAMstrandField intronMotif \
    --outFilterIntronMotifs RemoveNoncanonical \
    --outFilterMatchNminOverLread 0.66 \
    --outFilterScoreMinOverLread  0.66 \
    --outFilterMatchNmin  0 \
    --peOverlapNbasesMin  0 \
    --outReadsUnmapped Fastx

mv "$tempDir"/"$sid"_STAR* "$outputDir"/alignments
rm "$tempDir"/"$sid"_SQP_L001_R*_001.fastq.gz

mv "$outputDir"/alignment/"$sid"_STAR_*Unmapped* "$outputDir"/
cd "$outputDir"/

rename "_STAR_Unmapped.out.mate1" "_nohost_SQP_L001_R1_001.fastq" *
rename "_STAR_Unmapped.out.mate2" "_nohost_SQP_L001_R2_001.fastq" *

module load fastqc-0.11.9-gcc-12.1.0

fastqc "$outputDir"/alignment/"$sid"_nohost_SQP_L001_R*_001.fastq -t 4
mv "$outputDir"/alignment/"$sid"_nohost_SQP_L001_R*_001*fastqc* "$outputDir"/qc

gzip "$outputDir"/alignment/"$sid"_nohost_SQP_L001_R*_001.fastq
