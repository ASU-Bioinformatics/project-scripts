#!/bin/bash

#SBATCH -p gpu
#SBATCH -q wildfire
#SBATCH -t 0-00:15:00

#SBATCH -A kawoodbu
#SBATCH -o slurm.%j.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err               # STDERR (%j = JobId)
#SBATCH -C V100                       # bonito requires at least a volta level gpu card
#SBATCH ==gres=gpu:1

INPUT_FAST5_DIR="bonito-in"
OUTPUT_BASECALL_DIR="bonito-out"
OUTPUT_NAME_BASECALL="test.fasta"

#cd /data/gencore/analysis_projects/5609497_Dhenugen

cd /scratch/kawoodbu

module load cuda/11.3.0

source /data/biocore/programs/bonito/venv3.7/bin/activate

#module load tensorflow/1.8-agave-gpu
#module load tensorflow/1.12-py3
#source activate tf1.12-gpu

#python /data/biocore/programs/fast-bonito/basecaller.py --reads_directory $INPUT_FAST5_DIR \
#                                                        --output /$OUTPUT_BASECALL/${OUTPUT_NAME_BASECALL}

bonito basecaller dna_r10.4.1_e8.2_260bps_hac@v3.5.2 /data/gencore/analysis_projects/5609497_Dhenugen/fast5 > dhenugen-hacs.fastq

/data/biocore/programs/bonito/bonito/models/configs/
