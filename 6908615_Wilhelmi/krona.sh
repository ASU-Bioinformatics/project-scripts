#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 2-0:00
#SBATCH --mem=64G

module purge
METADATA="/data/gencore/analysis_projects/6908615_Wilhelmi_16s/metadata.txt"
FASTQ_DIR="/data/gencore/analysis_projects/6908615_Wilhelmi_16s/long-quality-trim/fastq"
OUTPUT_DIR="/data/gencore/analysis_projects/6908615_Wilhelmi_16s/qiime2"

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2023.7/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2023.7/gg_2022_10_backbone.v4.nb.qza"
#classifiers[18]="/data/biocore/qiime2_classifiers/qiime-2023.7/silva-138-99-nb-classifier.qza"
#classifiers[it]="/data/biocore/qiime2_classifiers/qiime-2022.2/unite-ver9-classifier-99-27-10-2022.qza"

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/qiime2023.7

cd "$OUTPUT_DIR"

for i in ${!classifiers[@]};
do

(echo "$i"

qiime krona collapse-and-plot \
--i-table "$OUTPUT_DIR"/table.qza \
--i-taxonomy "$OUTPUT_DIR"/"$i"-taxonomy.qza \
--o-krona-plot "$OUTPUT_DIR"/"$i"-krona.qzv
cd ..) &

done;