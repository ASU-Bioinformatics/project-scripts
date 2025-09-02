#!/bin/bash

##### check rRNA levels in RNAseq fastq output #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 0-12:00                         # estimated time needed
#SBATCH --mem=64G

module load mamba/latest

scriptsDir="/data/gencore/shared_scripts/RNAseq/old-scripts"
refSeq="/data/gencore/databases/reference_genomes/bovine/bos.taurus.rRNA.fna"

fastqDir="/data/gencore/analysis_projects/8424037_Molehin/fastq"
bbstatsDir="/data/gencore/analysis_projects/8424037_Molehin/bbstats"

mkdir -p "$bbstatsDir"

source activate /data/biocore/programs/conda-envs/bb-env

cd $fastqDir

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  echo $i
  bbduk.sh in="$i"_SRN_L001_R1_001.fastq.gz in2="$i"_SRN_L001_R2_001.fastq.gz \
           outm1="$bbstatsDir"/ribo-"$i"_SRN_L001_R1_001.fastq outm2="$bbstatsDir"/ribo-"$i"_SRN_L001_R2_001.fastq \
           outu1="$bbstatsDir"/nonribo-"$i"_SRN_L001_R1_001.fastq outu2="$bbstatsDir"/nonribo-"$i"_SRN_L001_R2_001.fastq \
           ref=$refSeq stats="$bbstatsDir"/"$i"_bbdukstats.txt
done;

source deactivate

cd $bbstatsDir
chmod -R g+w *
