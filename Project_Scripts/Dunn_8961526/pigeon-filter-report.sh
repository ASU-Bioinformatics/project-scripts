#!/bin/bash

##### pigeon filter and report #####
#### %j is job id

#SBATCH -p public         # sol
#SBATCH -q public         # sol
####SBATCH -p general        # phx
####SBATCH -q grp_kawoodbu   # phx
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-00:30
#SBATCH -c 1
#SBATCH --mem=16G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/ioseq-env

#pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.sorted.gff

#pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt
