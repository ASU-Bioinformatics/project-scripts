#!/bin/bash

##### qiime2 pt3: ancom #####

#SBATCH -p general
####SBATCH -q grp_kawoodbu
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 1-0:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=64G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7

fastqDir="/data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args/fastq"
qiimeDir="/data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2"
ancomDir="/data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2-ancombc"

mkdir -p $ancomDir

metadata="/data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt"

# sometimes these modules don't work if the --help call isn't run first
# it must be a bug
qiime dada2 --help
qiime composition --help

##### ANCOM with Taxa #####

### Classify, Visualize, and Collapse Table with 2024.10 at Class and Genus level ###
# can add different hierarchical levels!

## Formula and reference level need to be manually specified
# I change the file name to reflect these values
# especially when ancombc needs to be run for multiple categories

for i in 1 2 3 4 5 6 7;
do
  qiime composition ancombc \
    --i-table "$qiimeDir"/level"$i"_table-paired-paired.qza \
    --m-metadata-file "$metadata" \
    --p-formula Treatment \
    --p-reference-levels Treatment::control \
    --o-differentials "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza

  qiime composition da-barplot \
    --i-data "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza \
    --p-effect-size-label 'lfc' \
    --p-feature-id-label 'id' \
    --p-error-label 'se' \
    --p-significance-label 'q_val' \
    --p-significance-threshold 0.05 \
    --p-effect-size-threshold 0.0 \
    --p-level-delimiter ';' \
    --o-visualization "$ancomDir"/dabarplot-ancombc-treatment.control-level"$i"_table-paired-paired_qval.05.qzv

  qiime composition tabulate \
    --i-data "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza \
    --o-visualization "$ancomDir"/table-ancombc-treatment.control-level"$i"_table-paired-paired.qzv

  qiime tools export \
    --input-path "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza \
    --output-path "$ancomDir"/ancombc-treatment.control-level"$i"
done
