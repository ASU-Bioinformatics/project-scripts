#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 1-0:00                     # estimated time needed
#SBATCH -c 4
#SBATCH --mem=16G

module load bbmap-39.01-gcc-12.1.0

cd /data/gencore/analysis_projects/6078853_Otak_DNA/megahit-alignments/samples

for sid in $(find ./ -type f -name "*headed.sam" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
do
  pileup.sh in="$sid"_SQP_full.headed.sam \
          out="$sid"_covstats.txt \
          rpkm="$sid"_rpkm.txt
