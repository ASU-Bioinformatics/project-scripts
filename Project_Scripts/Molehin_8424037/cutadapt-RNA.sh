#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 2-0:00
#SBATCH --mem=64G

# this is the script I used for removing adapters and trimming for quality for Barrila RNA
# delete or comment out the merging code at the end if you're running it as an sbatch haha

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/cutadapt/

#mkdir -p /data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/cutadapt
cd /data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo

#gzip *.fastq

adapters="/data/gencore/databases/trimmomatic/PolyAndIllumina.fa"

#for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
#do
#  cutadapt -a file:"$adapters" -A file:"$adapters" \
#         -m 30 -q 20 \
#         -o "$i"-cut_SCT_L001_R1_001.fastq.gz \
#         -p "$i"-cut_SCT_L001_R2_001.fastq.gz \
#         "$i"_S*_L001_R1_001.fastq.gz "$i"_S*_L001_R2_001.fastq.gz
#done

#mv *cut* ./cutadapt

cd ./cutadapt

module load trimmomatic-0.39-gcc-12.1.0

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq)
do
  echo $i
  trimmomatic PE "$i"_SCT_L001_R1_001.fastq.gz "$i"_SCT_L001_R2_001.fastq.gz \
            "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
            "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
            CROP:150 SLIDINGWINDOW:4:15 MINLEN:36
done

mkdir -p ../cut-paired-filter0415
mkdir -p ../cut-unpaired-filter0415

mv *SQP* ../cut-paired-filter0415/
mv *SUN* ../cut-unpaired-filter0415/

cd ../
chmod -R g+w *
