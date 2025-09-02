#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8:00
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.UMD.out
#SBATCH -e slurm.%j.starGG.UMD.err


module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/bovine/bos-taurus-UMD3.1"
fasta="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna"
gff="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff"
gtf="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gtf"

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

chmod -R g+w *
#--sjdbGTFtagExonParentTranscript transcript \

### STAR Dual Genome Indexes ###
