#!/bin/bash

##### build bowtie2 index for ncbi viral database #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.E4.out              # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.E4.err              # STDERR (%j = JobId)
#SBATCH -c 16
#SBATCH -t 2-0:00
#SBATCH --mem=100G

module purge
module load mamba/latest

# if a bowtie index for the NCBI viral genome set already exists, this script should be skipped

##### Define Variables #####

dir="/data/gencore/databases/ncbi_viral"
sid="ncbi_viral"
assembly="$dir"/viral.1.1.genomic.fna

#### vironomy specific variables ####

dir="/data/gencore/analysis_projects/6078853_Otak_DNA/vibrant/VIBRANT_final.contigs/VIBRANT_phages_final.contigs"
sid="vibrant.predictions"
assembly="$dir/"final.contigs.phages_combined.fna
#### build bowtie2 index ####

source activate /data/biocore/programs/conda-envs/bowtie2-env

bowtie2-build --threads 64 "$assembly" "$dir"/"$sid"

source deactivate
