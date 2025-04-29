#!/bin/bash

##### fastqc file generation #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -t 2-0:00                    # estimated time needed
#SBATCH --mem=64G

module load trimmomatic-0.39-gcc-12.1.0
indir="/data/gencore/analysis_projects/6402395_Qian_16S/fastq/original"

#for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq);
#do
#  echo $i;
#  cat ./"$i"_S*_L001_R1_001.fastq.gz ./"$i"_S*_L002_R1_001.fastq.gz >> ./merged/"$i"_SME_L001_R1_001.fastq.gz
#  cat ./"$i"_S*_L001_R2_001.fastq.gz ./"$i"_S*_L002_R2_001.fastq.gz >> ./merged/"$i"_SME_L001_R2_001.fastq.gz
#done

#cd ./merged

cd $indir

for i in $(find ./ -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq);
do
  echo $i;
  trimmomatic PE -trimlog triminfo.txt \
              "$i"_S*_L00*_R1_001.fastq.gz "$i"_S*_L00*_R2_001.fastq.gz \
              "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
              "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
              ILLUMINACLIP:/data/gencore/databases/trimmomatic/TruSeq3-PE-2.fa:2:30:10 \
              LEADING:3 TRAILING:3 SLIDINGWINDOW:4:10 MINLEN:240
done
