#!/bin/bash

###SBATCH -p debug
###SBATCH -q wildfire
###SBATCH -t 0-00:15:00

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 7-0:00                     # estimated time needed
#SBATCH -c 4
#SBATCH --mem=32G

sid="$1"
assembly="$2"
fqDir="$3"
outDir="$4"

echo "bowtie alignment on $assembly for sample $sid"

module purge
module load mamba/latest

source activate /data/biocore/programs/conda-envs/bowtie2-env/
export PERL5LIB=/data/biocore/programs/conda-envs/bowtie2-env/bin/perl

cd "$outDir"

# use SME for untrimmed data, SQP for trimmed data
# also check lane ID. It's usually L001.

bowtie2 --verbose -p 8 -x "$assembly" \
        -1 "$fqDir"/"$sid"_SQP_L001_R1_001.fastq.gz \
        -2 "$fqDir"/"$sid"_SQP_L001_R2_001.fastq.gz \
        --very-sensitive -k 16 --np 1 --mp "1,1" --rdg "0,1" --rfg "0,1" \
        --score-min "L,0,-0.05" --no-unal -S "$outDir"/"$sid"_SQP_full.sam

#echo "full.sam file is complete"

module load samtools-1.16-gcc-11.2.0

tail -n+5 "$sid"_SQP_full.sam > "$sid"_SQP_full.headed.sam

echo "full.headed.sam file is complete"

samtools view -bS "$sid"_SQP_full.headed.sam > "$sid"_SQP_full.bam

echo "full.bam file is complete"

samtools sort "$sid"_SQP_full.bam -o "$sid"_SQP_full.sorted.bam

echo "full.sorted.bam file is complete"

samtools assembly "$sid"_SQP_full.sorted.bam

echo "full.sorted.bam is assembled"

samtools flagstat "$sid"_SQP_full.sorted.bam > "$sid"_flagstats.txt

# obtain coverage stats and RPKM values for all alignments

module load bbmap-39.01-gcc-12.1.0

pileup.sh in="$outDir"/"$sid"_SQP_full.sorted.bam \
          out="$outDir"/"$sid"_covstats.txt \
          rpkm="$outDir"/"$sid"_rpkm.txt
