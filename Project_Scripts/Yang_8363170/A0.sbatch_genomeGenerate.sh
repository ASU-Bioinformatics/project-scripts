#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.ecoli.out
#SBATCH -o slurm.%j.starGG.ecoli.err


module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2"
fasta="$refDir"/"CFT073.genomic.fna"
gtf="$refDir"/"CFT073.genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon CDS \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 3000000000000
