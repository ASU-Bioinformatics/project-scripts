#!/bin/bash

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%A.spades-reassemble.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.spades-reassemble.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 1-0:00
#SBATCH -c 4
#SBATCH --mem=128G

module load mamba/latest

python /data/biocore/programs/SPAdes/SPAdes-3.15.5-Linux/bin/spades.py --careful \
  -1 /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/reassembled_bins/reads_for_reassembly/bin.1.permissive_1.fastq \
  -2 /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/reassembled_bins/reads_for_reassembly/bin.1.permissive_2.fastq \
  --trusted-contigs /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/refined_bins/metawrap_70_5_bins/bin.1.fa \
  -o /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/spades_reassembled_bins/ \
  --memory 128 --threads 32
