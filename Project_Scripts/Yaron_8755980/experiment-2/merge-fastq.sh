#!/bin/bash

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%A.out
#SBATCH -e slurm.%A.err
#SBATCH -t 0-8:00
#SBATCH --mem=64G
#SBATCH -c 1

umask 0007

# experiment 2 merger
#dir1="/data/gencore/analysis_projects/8755980_Pallod/rerun-fastq/fastq-exp2"
#dir2="/data/gencore/analysis_projects/8755980_Pallod/experiment-2/unprocessed-fastq/fastq"
#outDir="/data/gencore/analysis_projects/8755980_Pallod/merged-fastq-exp2"

# experiment 1 merger
dir1="/data/gencore/analysis_projects/8755980_Pallod/rerun-fastq/fastq-exp1"
dir2="/data/gencore/sftp/j_yaron/8755980_Pallod_RNA-Experiment-1/1.unprocessed-fastq/fastq/"
outDir="/data/gencore/analysis_projects/8755980_Pallod/merged-fastq-exp1"

mkdir -p $outDir

for i in $(find $dir1 $dir2 -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq)
  do
    echo $i
    cat "$dir1"/"$i"_SRN_L00*_R1_001.fastq.gz "$dir2"/"$i"_SRN_L00*_R1_001.fastq.gz >> "$outDir"/"$i"_SME_L001_R1_001.fastq.gz
    cat "$dir1"/"$i"_SRN_L00*_R2_001.fastq.gz "$dir2"/"$i"_SRN_L00*_R2_001.fastq.gz >> "$outDir"/"$i"_SME_L001_R2_001.fastq.gz
done;

module load fastqc-0.12.1-gcc-11.2.0

cd $outDir
fastqc -t 256 *
mkdir ./fastq
mkdir ./qc
mv *fastq.gz ./fastq
mv *fastqc* ./qc

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/multiqc.v1.20/

multiqc $outDir/qc

mv multiqc* $outDir
