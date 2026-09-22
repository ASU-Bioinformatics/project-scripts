#!/bin/bash

##### pigeon classify #####
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

#pigeon classify \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.sorted.gff \
#  /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.sorted.gtf \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  --fl /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.flnc_count.txt

pigeon classify \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.sorted.gff \
  /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.sorted.gtf \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  --fl /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.flnc_count.txt

pigeon classify \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.sorted.gff \
  /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.sorted.gtf \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  --fl /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.flnc_count.txt

pigeon classify \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.sorted.gff \
  /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.sorted.gtf \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  --fl /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.flnc_count.txt
