#!/bin/bash

##### pbmm2 test #####
#### %j is job id

#SBATCH -p general
#SBATCH -q grp_kawoodbu
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-04:00
#SBATCH -c 6
#SBATCH --mem=64G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/ioseq-env

#isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc03_5p--IsoSeqX_3p.bam \
#              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
#              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p.polya.flnc.bam

#isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc04_5p--IsoSeqX_3p.bam \
#              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
#              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p.polya.flnc.bam

#isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc11_5p--IsoSeqX_3p.bam \
#              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
#              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p.polya.flnc.bam

isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc12_5p--IsoSeqX_3p.bam \
              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p.polya.flnc.bam
