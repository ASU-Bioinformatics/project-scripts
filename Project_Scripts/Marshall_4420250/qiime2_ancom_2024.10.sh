#!/bin/bash

##### qiime2 pt3: ancom #####

#SBATCH -p general
####SBATCH -q grp_kawoodbu
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 1-0:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=64G

#### the ancom module appears to have been updated a lot since this version of qiime2, so I'm also wanting to test it with the newest verion
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2024.10-try2

DIRECTORY="/data/gencore/analysis_projects/4420250_Marshall_Publication"
FASTQ_DIR="$DIRECTORY"/fastq
OUTPUT_DIR="$DIRECTORY"/qiime2

METADATA="$DIRECTORY"/metadata.txt

qiime dada2 --help
qiime composition --help

#qiime tools import \
#  --type 'SampleData[PairedEndSequencesWithQuality]' \
#  --input-path "$FASTQ_DIR" \
#  --input-format CasavaOneEightSingleLanePerSampleDirFmt \
#  --output-path "$OUTPUT_DIR"/demux-paired-end-2024.10.qza

# DADA2: trim-length depends on quality visualization. this quality control process will additionally filter any phiX reads (commonly present in marker gene Illumina sequence data) that are identified in the sequencing data, and will filter chimeric sequences.
#qiime dada2 denoise-single \
#  --i-demultiplexed-seqs "$OUTPUT_DIR"/demux-paired-end-2024.10.qza \
#  --p-trim-left 0 \
#  --p-trunc-len 240 \
#  --p-n-threads 24 \
#  --o-representative-sequences "$OUTPUT_DIR"/rep-seqs-2024.10-sol.qza \
#  --o-table "$OUTPUT_DIR"/table-2024.10-sol.qza \
#  --o-denoising-stats "$OUTPUT_DIR"/stats-dada2-2024.10-sol.qza

#qiime composition add-pseudocount \
#  --i-table "$OUTPUT_DIR"/table-2024.10-sol.qza \
#  --o-composition-table "$OUTPUT_DIR"/comp-table-2024.10-sol.qza

# ancombc analysis (formula-based)
#qiime composition ancombc \
#  --i-table "$OUTPUT_DIR"/table-2024.10-sol.qza \
#  --m-metadata-file "$METADATA" \
#  --p-formula Replicate \
#  --o-differentials "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol.qza

#qiime composition ancombc \
#  --i-table "$OUTPUT_DIR"/table-2024.10-sol.qza \
#  --m-metadata-file "$METADATA" \
#  --p-formula 'Replicate + Location + treatment' \
#  --o-differentials "$OUTPUT_DIR"/ancombc-all-2024.10-sol.qza

# ancombc tabular output visualization
qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol.qza \
  --o-visualization "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol.qzv

#qiime composition tabulate \
#  --i-data "$OUTPUT_DIR"/ancombc-all-2024.10-sol.qza \
#  --o-visualization "$OUTPUT_DIR"/ancombc-all-2024.10-sol.qzv

# ancombc barplot output visualization
#qiime composition da-barplot \
#  --i-data "$OUTPUT_DIR"/ancombc-all-2024.10-sol.qza \
#  --p-effect-size-label 'lfc' \
#  --p-feature-id-label 'id' \
#  --p-error-label 'se' \
#  --p-significance-label 'q_val' \
#  --p-significance-threshold 0.05 \
#  --p-effect-size-threshold 0.0 \
#  --o-visualization "$OUTPUT_DIR"/ancombc-all-2024.10-sol-qval.05.qza

#qiime composition da-barplot \
#  --i-data "$OUTPUT_DIR"/ancombc-all-2024.10-sol.qza \
#  --p-effect-size-label 'lfc' \
#  --p-feature-id-label 'id' \
#  --p-error-label 'se' \
#  --p-significance-label 'q_val' \
#  --p-significance-threshold 1.0 \
#  --p-effect-size-threshold 0.0 \
#  --o-visualization "$OUTPUT_DIR"/ancombc-all-2024.10-sol-qval1.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --o-visualization "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol-qval.05.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 1.0 \
  --p-effect-size-threshold 0.0 \
  --o-visualization "$OUTPUT_DIR"/ancombc-replicate-2024.10-sol-qval1.qzv

# ancom analysis (without bias correction, with visualization)
#qiime composition ancom \
#  --i-table "$OUTPUT_DIR"/table-2024.10-sol.qza \
#  --m-metadata-file "$METADATA" \
#  --m-metadata-column 'Replicate' \
#  --p-transform-function 'clr' \
#  --p-difference-function 'f_statistic' \
#  --o-visualization "$OUTPUT_DIR"/ancom-replicate-2024.10-sol.qzv

##### ANCOM with Taxa #####
#classifier="/data/biocore/qiime2_classifiers/qiime-2024.5/silva-138-99-nb-classifier.qza"

#qiime feature-classifier classify-sklearn \
#  --i-classifier $classifier \
#  --i-reads rep-seqs-2024.10-sol.qza \
#  --o-classification sv-taxonomy-2024.10-sol.qza

#qiime metadata tabulate \
#  --m-input-file sv-taxonomy-2024.10-sol.qza \
#  --o-visualization sv-taxonomy-2024.10-sol.qzv

#qiime taxa barplot \
#  --i-table table-2024.10-sol.qza \
#  --i-taxonomy sv-taxonomy-2024.10-sol.qza \
#  --m-metadata-file "$METADATA" \
#  --o-visualization taxa-bar-plots-2024.10-sol.qzv

#qiime taxa collapse \
#  --i-table table-2024.10-sol.qza \
#  --i-taxonomy sv-taxonomy-2024.10-sol.qza \
#  --o-collapsed-table level3_table-2024.10-sol.qza \
#  --p-level 3

#qiime feature-table relative-frequency \
#  --i-table level3_table-2024.10-sol.qza \
#  --o-relative-frequency-table rel-level3_table-2024.10-sol.qza

#qiime taxa collapse \
#  --i-table table-2024.10-sol.qza \
#  --i-taxonomy sv-taxonomy-2024.10-sol.qza \
#  --o-collapsed-table level6_table-2024.10-sol.qza \
#  --p-level 6

#qiime feature-table relative-frequency \
#--i-table level6_table-2024.10-sol.qza \
#--o-relative-frequency-table rel-level6_table-2024.10-sol.qza

qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/level3_table-2024.10-sol.qza \
  --m-metadata-file "$METADATA" \
  --p-formula Replicate \
  --p-reference-levels Replicate::snsp \
  --o-differentials "$OUTPUT_DIR"/ancombc-replicate.snsp-level3-2024.10-sol.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.snsp-level3-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-replicate.snsp-level3-2024.10-sol-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.snsp-level3-2024.10-sol.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-replicate.snsp-level3-2024.10-sol.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-replicate.snsp-level3-2024.10-sol.qza \
  --output-path "$OUTPUT_DIR"/ancombc-snsp-level3


qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/level3_table-2024.10-sol.qza \
  --m-metadata-file "$METADATA" \
  --p-formula Replicate \
  --p-reference-levels Replicate::snpe \
  --o-differentials "$OUTPUT_DIR"/ancombc-replicate.snpe-level3-2024.10-sol.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.snpe-level3-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-replicate.snpe-level3-2024.10-sol-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.snpe-level3-2024.10-sol.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-replicate.snpe-level3-2024.10-sol.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-replicate.snpe-level3-2024.10-sol.qza \
  --output-path "$OUTPUT_DIR"/snpe

## reference level SUSP
qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/level3_table-2024.10-sol.qza \
  --m-metadata-file "$METADATA" \
  --p-formula Replicate \
  --p-reference-levels Replicate::susp \
  --o-differentials "$OUTPUT_DIR"/ancombc-replicate.susp-level3-2024.10-sol.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.susp-level3-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-replicate.susp-level3-2024.10-sol-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.susp-level3-2024.10-sol.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-replicate.susp-level3-2024.10-sol.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-replicate.susp-level3-2024.10-sol.qza \
  --output-path "$OUTPUT_DIR"/ancombc-susp-level3

## reference level SPE
qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/level3_table-2024.10-sol.qza \
  --m-metadata-file "$METADATA" \
  --p-formula Replicate \
  --p-reference-levels Replicate::spe \
  --o-differentials "$OUTPUT_DIR"/ancombc-replicate.spe-level3-2024.10-sol.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.spe-level3-2024.10-sol.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-replicate.spe-level3-2024.10-sol-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-replicate.spe-level3-2024.10-sol.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-replicate.spe-level3-2024.10-sol.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-replicate.spe-level3-2024.10-sol.qza \
  --output-path "$OUTPUT_DIR"/ancombc-spe-level3
