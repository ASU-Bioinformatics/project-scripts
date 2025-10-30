#!/bin/bash

##### lima demultiplexing of Kinnex bam files #####

#SBATCH -o slurm.%j.lima.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.lima.err               # STDERR (%j = JobId)
#SBATCH -p general
#SBATCH -q public
#SBATCH -t 1-0:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/lima-env

cd /data/gencore/analysis_projects/PacBio-NotKinnex-20250206/try3

# because lima is erroring out (segmentation fault, core dumped) when using the bam as an input file
# running bam2fastq first lets us run lima demultiplexing with fastq input and avoid that error
bam2fastq -o m84132-combined m84132_250124_184858_s3.hifi_reads.bam

# this primer file contains the long barcode sequences (the unique index plus the additional conserved adjacent region)
# using only the unique index failed to capture the majority of the barcode pairs
# the qc filter settings used are less stringent than the hifi preset, because more than 95% of the reads were lost with the presets
lima m84132-combined.fastq.gz \
     /data/gencore/analysis_projects/completed_projects/Kinnex-Project/16s_primers.fasta \
     m84132-combined-demux-from-fastq.fastq --split \
     --peek-guess --ccs --min-score 50 --min-end-score 30 \
     --min-ref-span 0.5 --different --min-scoring-regions 2

# rename the fastq files to match the sample names from the metadata

mkdir renamed-fastq
while IFS=, read kinnex sample; do
  cp $kinnex ./renamed-fastq/$sample
done < kinnex-names-try3.csv
