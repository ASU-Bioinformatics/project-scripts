#!/bin/bash

#### fastqc file generation #####

#SBATCH -p general
#SBATCH -q public #sol
#SBATCH -t 1-0:00               # estimated time needed
#SBATCH --mem=256G
#SBATCH -c 2

umask 0007
module purge
module load fastqc-0.12.1-gcc-11.2.0 #sol

##### go to directory where fastq files are located #####

#FDIR="/data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/nonribo"
#FDIR="/data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/cutadapt"
FDIR="/data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/cut-paired-filter0415"

cd "$FDIR"
mkdir fastq
mkdir qc

fastqc -t 256 *.fastq.gz
mv *fastqc* qc
mv *fastq.gz fastq

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/multiqc.v1.20/
cd "$FDIR"
multiqc "$FDIR/qc"

source deactivate

chmod -R g+w *
