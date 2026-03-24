#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p public
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8:00
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.mouse-norRNA.out
#SBATCH -e slurm.%j.starGG.mouse-norRNA.err


module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/GRCm39-no-rRNA"
fasta="$refDir"/"Mus_musculus.GRCm39.dna.primary_assembly.fa"
gtf="$refDir"/"Mus_musculus.GRCm39.114.no-rRNA.gtf"

cd "$refDir"

#gffread "$gff" -T -o "$gtf"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon exon \
  --limitGenomeGenerateRAM 3000000000000

chmod -R g+w *
#--sjdbGTFtagExonParentTranscript transcript \

### STAR Dual Genome Indexes ###
