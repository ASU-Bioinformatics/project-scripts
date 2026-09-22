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

# using the Revio output files:

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc03_5p--IsoSeqX_3p.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc03_5p--IsoSeqX_3p.bam

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc04_5p--IsoSeqX_3p.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc04_5p--IsoSeqX_3p.bam

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc11_5p--IsoSeqX_3p.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc11_5p--IsoSeqX_3p.bam

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/ReUploadBC12/IsoSeqX_bc12_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc12_5p--IsoSeqX_3p.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/ReUploadBC12/IsoSeqX_bc12_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc12_5p--IsoSeqX_3p.bam

# using the clustered bam files:

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.clustered.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.aligned.bam

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.clustered.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.aligned.bam

#pbmm2 align --preset ISOSEQ \
#  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.clustered.bam \
#  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
#  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.aligned.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.clustered.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.aligned.bam
