#!/bin/bash

##### iso-collapse #####
#### %j is job id

#SBATCH -p general
#SBATCH -q grp_kawoodbu
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-04:00
#SBATCH -c 2
#SBATCH --mem=32G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/ioseq-env

#isoseq collapse --do-not-collapse-extra-5exons \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.aligned.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.polya.flnc.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.gff

#isoseq collapse --do-not-collapse-extra-5exons \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.aligned.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.polya.flnc.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.gff

#isoseq collapse --do-not-collapse-extra-5exons \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.aligned.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.polya.flnc.bam \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.gff

isoseq collapse --do-not-collapse-extra-5exons \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.aligned.bam \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.polya.flnc.bam \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.gff
