#!/bin/bash

##### align metatranscriptomic reads to paired metagenomic predicted genes #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.E5.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.E5.err               # STDERR (%j = JobId)
#SBATCH -t 0-00:15

module purge
module load mamba/latest

source activate /data/biocore/programs/conda-envs/bowtie2-env/
export PERL5LIB=/data/biocore/programs/conda-envs/bowtie2-env/bin/perl

assemblyDir="/data/gencore/analysis_projects/6078853_Otak_DNA/vibrant/VIBRANT_final.contigs/VIBRANT_phages_final.contigs"
fastqDir="/data/gencore/analysis_projects/6078853_Otak_DNA/fastq/adapter-trimmed/fastq"
alignmentDir="/data/gencore/analysis_projects/6078853_Otak_DNA/vibrant/vibrant-alignments"
asmPrefix="vibrant.predictions"
use="DIRECTORY"
list=""

mkdir -p "$alignmentDir"

if [ "$use" == "DIRECTORY" ];
then

### run E3.bowtie-align for all samples in a directory ###

for i in $(find "$fastqDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
do
  echo "$i"
  sbatch --job-name "$i".E5.bowtie-align \
       /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_E/E5.bowtie-align.sh \
       "$i" "$assemblyDir" "$fastqDir" "$alignmentDir" "$asmPrefix"
done;

elif [ "$use" == "LIST" ];
then

### run D3.bowtie-align for a subset of samples in a directory ###

for i in $list
do
  echo "$i"
  sbatch --job-name "$i".E5.bowtie-align \
       /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_E/E5.bowtie-align.sh \
       "$i" "$asmPrefix" "$contigDir" "$fqDir" "$alignDir"
done;

fi
