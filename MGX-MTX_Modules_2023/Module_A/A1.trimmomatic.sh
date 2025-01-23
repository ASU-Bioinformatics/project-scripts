#!/bin/bash

#### Run on Sol ####

#SBATCH -p general
#SBATCH -q public
#SBATCH -t 2-0:00
#SBATCH --mem=64G

#### Load Modules ####

module load trimmomatic-0.39-gcc-12.1.0
module load fastqc-0.11.9-gcc-12.1.0

#### Define Variables ####

projectDir=$1
cd "$projectDir"

#### Trim and Filter Reads ####
for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq);
do
  echo $i;
  trimmomatic PE "$i"_SRN_L00*_R1_001.fastq.gz "$i"_SRN*_L00*_R2_001.fastq.gz \
              "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
              "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
              ILLUMINACLIP:/data/gencore/databases/trimmomatic/TruSeq3-PE-2.fa:2:30:10 \
              CROP:151 SLIDINGWINDOW:4:15 MINLEN:100
done

#### Organize Files ####

mkdir unpaired-lqt
mkdir adapter-trimmed
mkdir original
mv *SQP* adapter-trimmed/
mv *SUN* unpaired-lqt/
mv *fastq.gz original

#### Fastqc Trimmed Data ####

cd adapter-trimmed/
fastqc -t 16 *
mkdir fastq
mkdir qc
mv *fastq.gz fastq/
mv *fastqc* qc/

# typically will run MultiQC manually since it is quick and can get confusing if any errors occurred upstream
