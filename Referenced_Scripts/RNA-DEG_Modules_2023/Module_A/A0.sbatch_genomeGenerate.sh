#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public
#SBATCH -t 0-8
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.out
#SBATCH -o slurm.%j.starGG.err


module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/impatiens_spp/impatiens-glandulifera-simplified-labels"
fasta="$refDir"/"glandulifera-9.labeled.fasta"
gff="$refDir"/"glandulifera-9.labeled.gff"
gtf="$refDir"/"glandulifera-9.labeled.gtf"

cd "$refDir"

gffread "$gff" -T -o "$gtf"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon exon \
  --limitGenomeGenerateRAM 3000000000000


#--sjdbGTFtagExonParentTranscript transcript \

### STAR Dual Genome Indexes ###
