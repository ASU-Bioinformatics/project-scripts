#!/bin/bash

##### pigeon prepare #####
#### %j is job id

#SBATCH -p public         # sol
#SBATCH -q public         # sol
####SBATCH -p general        # phx
####SBATCH -q grp_kawoodbu   # phx
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-03:30
#SBATCH -c 1
#SBATCH --mem=16G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/ioseq-env

#pigeon prepare \
  #/data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.gff \
  #/data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.gff \
  #/data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.gff \
  #/data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.gff \


pigeon prepare /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.gtf
