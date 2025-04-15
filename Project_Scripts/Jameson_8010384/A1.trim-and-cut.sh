#!/bin/bash

#### Run on Sol ####

#SBATCH -p general
#SBATCH -q public
#SBATCH -t 2-0:00
#SBATCH --mem=64G
#SBATCH -o slurm.%j.A1.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.A1.err               # STDERR (%j = JobId)

#### Load Modules and Cutadapt Environment ####

module load trimmomatic-0.39-gcc-12.1.0
module load fastqc-0.11.9-gcc-12.1.0
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/cutadapt/

#### Define Variables ####

projectDir="/data/gencore/analysis_projects/8010384_coral-MGX/fastq"
cd "$projectDir"

adapters="/data/gencore/databases/trimmomatic/TruSeq3-PE-2.fa"

#### Remove Adapter Sequence ####
for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq)
do
  cutadapt -a file:"$adapters" -A file:"$adapters" \
         -m 30 -q 20 \
         -o "$i"-cut_SCT_L001_R1_001.fastq.gz \
         -p "$i"-cut_SCT_L001_R2_001.fastq.gz \
         "$i"_S*_L001_R1_001.fastq.gz "$i"_S*_L001_R2_001.fastq.gz
done

#### Trim and Filter Reads ####
for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq);
do
  echo $i;
  trimmomatic PE "$i"_SCT_L00*_R1_001.fastq.gz "$i"_SCT*_L00*_R2_001.fastq.gz \
              "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
              "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
              CROP:151 SLIDINGWINDOW:4:15 MINLEN:100
done

#### Organize Files ####

mkdir -p cutadapt-only
mkdir -p cut-unpaired0415
mkdir -p cut-paired0415
mkdir -p original
mv *SCT* cutadapt-only/
mv *SQP* cut-paired0415/
mv *SUN* cut-unpaired0415/
mv *fastq.gz original

#### Fastqc Trimmed Data ####

cd cut-paired0415/
fastqc -t 32 *
mkdir fastq
mkdir qc
mv *fastq.gz fastq/
mv *fastqc* qc/

# typically will run MultiQC manually since it is quick and can get confusing if any errors occurred upstream
