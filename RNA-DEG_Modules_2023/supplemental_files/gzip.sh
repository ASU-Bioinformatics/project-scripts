#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.gzip.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.gzip.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 2-0:00
#SBATCH -c 6
#SBATCH --mem=128G

for i in $(find ./ -type f -name "*.fastq" | while read F; do basename $F; done)
do
  gzip $i
done
