#!/bin/bash

##### classify and build phylogeny for viral contigs from vibrant #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH -o slurm.%j.E2.vironomy.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.E2.vironomy.err               # STDERR (%j = JobId)
#SBATCH -t 2-0:00
#SBATCH -c 16
#SBATCH --mem=256G

#### set up environment ####

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/vironomy-env

#### download databases (only do once) ####

#vironomy download_db -p /data/gencore/databases/vironomy

#### vironomy end-to-end (both classify and phylogeny)

contigs="/data/gencore/sftp/otakuye_conroyben/6078853_Metagenomics/10.viral-analysis/VIBRANT_final.contigs/final.contigs.phages_combined.fna"
dbDir="/data/gencore/databases/vironomy/vironomy_database"
outDir="/data/gencore/analysis_projects/6078853_Otak_DNA/vironomy-prodigal"

mkdir -p $outDir
cp $contigs $outDir
cd $outDir

sed "/^>/ s/ /_/g" < final.contigs.phages_combined.fna > newhead.final.contigs.phages_combined.fna

vironomy end_to_end -i newhead.final.contigs.phages_combined.fna \
                    -d /data/gencore/databases/vironomy/vironomy_database \
                    -p 16

vironomy phylogeny -c /data/gencore/analysis_projects/6078853_Otak_DNA/vironomy/vironomy_output \
                   -d /data/gencore/databases/vironomy/vironomy_database \
                   -o /data/gencore/analysis_projects/6078853_Otak_DNA/vironomy/phylogeny-rerun

# use ViPTree to manually construct phylogenies that incorporate viral databases with identified viruses
# separate contigs by phylum (or lower taxon if necessary) based on vironomy classification
# use seqtk subseq (biocore-rna environment) to extract contigs into subsetted fasta files
