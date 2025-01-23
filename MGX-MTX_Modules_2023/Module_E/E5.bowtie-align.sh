#!/bin/bash

###SBATCH -p debug
###SBATCH -q wildfire
###SBATCH -t 0-00:15:00

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 2-0:00                     # estimated time needed
#SBATCH -c 4
#SBATCH --mem=16G

#### read in variables ####
sid=$1
assemblyDir=$2
fastqDir=$3
alignmentDir=$4
asmPrefix=$5
echo "bowtie alignment on $asmPrefix fasta assembly for sample $sid"

module purge
module load mamba/latest

source activate /data/biocore/programs/conda-envs/bowtie2-env/
export PERL5LIB=/data/biocore/programs/conda-envs/bowtie2-env/bin/perl

cd "$alignmentDir"

bowtie2 --verbose -p 8 -x "$assemblyDir"/"$asmPrefix" \
        -1 "$fastqDir"/"$sid"_SQP_L001_R1_001.fastq.gz \
        -2 "$fastqDir"/"$sid"_SQP_L001_R2_001.fastq.gz \
        --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "L,0,-0.05" --no-unal -S "$alignmentDir"/"$sid"_SQP_full.sam

tail -n+5 "$alignmentDir"/"$sid"_SQP_full.sam > "$alignmentDir"/"$sid"_SQP_full.headed.sam

samtools view -bS "$alignmentDir"/"$sid"_SQP_full.headed.sam > "$alignmentDir"/"$sid"_SQP_full.bam

samtools sort "$alignmentDir"/"$sid"_SQP_full.bam -o "$alignmentDir"/"$sid"_SQP_full.sorted.bam

samtools index "$alignmentDir"/"$sid"_SQP_full.sorted.bam

# obtain coverage stats and RPKM values for all alignments

source activate /data/biocore/programs/conda-envs/bb-env

pileup.sh in="$alignmentDir"/"$sid"_SQP_full.headed.sam \
          out="$alignmentDir"/"$sid"_covstats.txt \
          rpkm="$alignmentDir"/"$sid"_rpkm.txt
