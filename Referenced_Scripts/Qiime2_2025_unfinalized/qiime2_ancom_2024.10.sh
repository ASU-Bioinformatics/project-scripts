#!/bin/bash

##### qiime2 pt3: ancom #####

#SBATCH -p general
####SBATCH -q grp_kawoodbu
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 1-0:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=64G

# the formulas and control levels need to be specified manually at this point
# my ultimate goal is to incorporate them into a parameter-style input
# however, all taxonomic levels will be evaluated for the input formula
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2024.10-try2

metadata="$1"
fastqDir="$2"
outDir="$3"

##### ANCOM with Taxa #####
classifier="/data/biocore/qiime2_classifiers/qiime-2024.5/silva-138-99-nb-classifier.qza"

cd "$outDir"

for i in 1 2 3 4 5 6 7;
do
  echo $i
  qiime taxa collapse \
    --i-table "$outDir"/table.qza \
    --i-taxonomy "$outDir"/sv-taxonomy.qza \
    --o-collapsed-table "$outDir"/level"$i"_table.qza \
    --p-level $i

  qiime feature-table relative-frequency \
    --i-table "$outDir"/level"$i"_table.qza \
    --o-relative-frequency-table "$outDir"/rel-level"$i"_table.qza

  qiime composition ancombc \
    --i-table "$outDir"/level"$i"_table.qza \
    --m-metadata-file "$METADATA" \
    --p-formula Replicate \
    --p-reference-levels Replicate::snsp \
    --o-differentials "$outDir"/ancombc-replicate.snsp-level"$i".qza

  qiime composition da-barplot \
    --i-data "$outDir"/ancombc-replicate.snsp-level3-2024.10-sol.qza \
    --p-effect-size-label 'lfc' \
    --p-feature-id-label 'id' \
    --p-error-label 'se' \
    --p-significance-label 'q_val' \
    --p-significance-threshold 0.05 \
    --p-effect-size-threshold 0.0 \
    --p-level-delimiter ';' \
    --o-visualization "$outDir"/dabarplot-ancombc-replicate.snsp-level"$i"-qval.05.qza

  qiime composition tabulate \
    --i-data "$outDir"/ancombc-replicate.snsp-level"$i".qza \
    --o-visualization "$outDir"/table-ancombc-replicate.snsp-level"$i".qzv

  qiime tools export \
    --input-path "$outDir"/ancombc-replicate.snsp-level3"$i".qza \
    --output-path "$outDir"/ancombc-snsp-level"$i"

done
