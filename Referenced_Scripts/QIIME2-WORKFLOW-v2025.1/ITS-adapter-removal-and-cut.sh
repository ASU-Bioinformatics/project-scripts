#!/bin/bash

#SBATCH -p public
#SBATCH -q public                     # sol QOS
####SBATCH -q grp_kawoodbu           # phx QOS
#SBATCH -o slurm.%A.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-12:00
#SBATCH --mem=32G

# this is the script is used for removing adapters and trimming for quality, specifically for 18S reads (so the reads are limited to a minimum of 30bp and a maximum of 150bp)
# the sbatch parameters are set for sol by default

#the script is designed to take any input fastq.gz files where the sample ID does not contain an underscore and is separated from the remainder of
# the file name by an underscore, as long as a R1 or R2 can be found in between the sample ID and extension.
# the output files, however, will be formatted with the traditional Illumina/Casava naming system.
umask 007

inputDir=$1
trimDir=$2
cutDir=$3

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/cutadapt/

mkdir -p $trimDir
mkdir -p $cutDir

cd $inputDir

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq )
do
  ( echo "$i"
    cutadapt \
      --nextseq-trim=12 \
      -g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGCTCTTGGTCATTTAGAGGAAGTAA;anywhere;min_overlap=6" \
      -G "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGATGCTGCGTTCTTCATCGATGC;anywhere;min_overlap=6" \
      -A "CTGTCTCTTATA;anywhere;min_overlap=10" \
      -a "CTGTCTCTTATA;anywhere;min_overlap=10" \
      -n 2 \
      -o "$trimDir"/"$i"_SC1_L001_R1_001.fastq.gz \
      -p "$trimDir"/"$i"_SC1_L001_R2_001.fastq.gz \
      "$i"_S*R1*.fastq.gz "$i"_S*R2*.fastq.gz

    cutadapt \
      -m 30 -l 250 \
      -o "$cutDir"/"$i"_SCT_L001_R1_001.fastq.gz \
      -p "$cutDir"/"$i"_SCT_L001_R2_001.fastq.gz \
      "$trimDir"/"$i"_SC1_L001_R1_001.fastq.gz "$trimDir"/"$i"_SC1_L001_R1_001.fastq.gz
  ) &
done;
wait

source deactivate
source activate /data/biocore/programs/mamba-envs/multiqc.v1.20

module load fastqc-0.12.1-gcc-11.2.0

cd $cutDir

fastqc -t 256 *
mkdir fastq
mkdir fastqc

mv *.fastq.gz fastq
mv *fastqc* fastqc

cd fastqc
multiqc *
mv multiqc* ../

cd $trimDir

fastqc -t 256 *
mkdir fastq
mkdir fastqc

mv *.fastq.gz fastq
mv *fastqc* fastqc

cd fastqc
multiqc *
mv multiqc* ../
