#!/bin/bash

##### qiime2 pt1: dada2 #####

#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 0-4:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=32G

# call this script as: "sbatch qiime2_pt1.sh metadata fastqDir outDir"
# there are enough optional settings for dada2 that I may up this control to a parameter style script...

module purge
metadata="$1"
fastqDir="$2"
outDir="$3"

mkdir -p "$outDir"

source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2024.10-try2

cd "$projectDir"

# Import the fastq data into the Qiime2 interface
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$fastqDir" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "$outDir"/demux-paired-end.qza

# Generate a summary of the demultiplexing results
qiime demux summarize \
  --i-data "$outDir"/demux-paired-end.qza \
  --o-visualization "$outDir"/demux.qzv

# DADA2: trim-length depends on quality visualization.
# This quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data)
# and will filter chimeric sequences.
# The paired-read input option is preferred with paired-end data;
# however, if dada2 is crashing or not letting reasonable amounts of data through its filters,
# single end settings can be used to avoid issues caused by insufficient overlap length.

qiime dada2 denoise-single \
  --i-demultiplexed-seqs "$outDir"/demux-paired-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 240 \
  --o-representative-sequences "$outDir"/rep-seqs.qza \
  --o-table "$outDir"/table.qza \
  --o-denoising-stats "$outDir"/stats-dada2.qza

# visualize metadata stats from denoising
qiime metadata tabulate \
  --m-input-file "$outDir"/stats-dada2.qza \
  --o-visualization "$outDir"/stats-dada2.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table "$outDir"/table.qza \
  --o-visualization "$outDir"/table.qzv \
  --m-sample-metadata-file "$metadata"

qiime feature-table tabulate-seqs \
  --i-data "$outDir"/rep-seqs.qza \
  --o-visualization "$outDir"/rep-seqs.qzv
