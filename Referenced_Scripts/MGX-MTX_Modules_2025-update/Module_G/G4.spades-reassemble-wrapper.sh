#!/bin/bash

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%A.spades-reassemble.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.spades-reassemble.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 1-0:00
#SBATCH -c 4
#SBATCH --mem=128G

module load mamba/latest

inputDir="/data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/refined_bins/metawrap_70_5_bins"
samplenum=$(ls $inputDir | wc -l)
script="/data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_G/G4A.spades-reassemble.sh"

sbatch --array=0-$((samplenum-1)) --parsable --export=ALL "$script"
