#!/bin/bash

##### qiime2 pt1: dada2 #####

#SBATCH -A kawoodbu
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH --mail-type=ALL                   # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=kawoodbu@asu.edu    # send-to address
#SBATCH -t 0-4:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=32G

module purge

DIRECTORY="/data/gencore/analysis_projects/5327771_Ball_16S"
FASTQ_DIR="$DIRECTORY"/fastq
OUTPUT_DIR="$DIRECTORY"/qiime-results
METADATA="$DIRECTORY"/scripts/"metadata.tsv"

source activate /data/biocore/programs/conda-envs/qiime2-2022.2

cd "$DIRECTORY"

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "$FASTQ_DIR" \
  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
  --output-path "$OUTPUT_DIR"/demux-paired-end.qza

# Generate a summary of the demultiplexing results
qiime demux summarize \
  --i-data "$OUTPUT_DIR"/demux-paired-end.qza \
  --o-visualization "$OUTPUT_DIR"/demux.qzv

# DADA2: trim-length depends on quality visualization. this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "$OUTPUT_DIR"/demux-paired-end.qza \
  --p-trim-left 0 \
  --p-trunc-len 240 \
  --o-representative-sequences "$OUTPUT_DIR"/rep-seqs.qza \
  --o-table "$OUTPUT_DIR"/table.qza \
  --o-denoising-stats "$OUTPUT_DIR"/stats-dada2.qza

# visualize metadata stats from denoising
qiime metadata tabulate \
  --m-input-file "$OUTPUT_DIR"/stats-dada2.qza \
  --o-visualization "$OUTPUT_DIR"/stats-dada2.qzv

qiime feature-table summarize \
  --i-table "$OUTPUT_DIR"/table.qza \
  --o-visualization "$OUTPUT_DIR"/table-quick-check.qzv

# FeatureTable and FeatureData summaries
qiime feature-table summarize \
  --i-table "$OUTPUT_DIR"/table.qza \
  --o-visualization "$OUTPUT_DIR"/table.qzv \
  --m-sample-metadata-file "$METADATA"

qiime tools export --input-path "$OUTPUT_DIR"/table.qza --output-path ./asv-table
cd ./asv-table
biom convert --to-tsv -i feature-table.biom -o asv-table.tsv
cd ..

#edit asv-table.tsv in R with microdecon

#biom convert -i decon-asvs.tsv -o decon-asvs.biom --table-type="OTU table" --to-json

#qiime tools import \
#  --input-path decon-asvs.biom \
#  --type 'FeatureTable[Frequency]' \
#  --input-format BIOMV100Format \
#  --output-path "$OUTPUT_DIR"/decon-table.qza

#qiime feature-table summarize \
#  --i-table "$OUTPUT_DIR"/decon-table.qza \
#  --o-visualization "$OUTPUT_DIR"/decon-table.qzv \
#  --m-sample-metadata-file "$METADATA"

qiime feature-table tabulate-seqs \
  --i-data "$OUTPUT_DIR"/rep-seqs.qza \
  --o-visualization "$OUTPUT_DIR"/rep-seqs.qzv
