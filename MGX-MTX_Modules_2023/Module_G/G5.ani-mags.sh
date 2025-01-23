#!/bin/bash

##### ani calculations for Impatiens dataset #####

#SBATCH -t 0-4:00                    # estimated time needed
#SBATCH --mem=72G
#SBATCH --cpus-per-task=16

module purge

source activate /data/biocore/programs/conda-envs/pyani-env/

module load mummer/4.0.0
module load blast/2.11.0

cd /data/gencore/analysis_projects/6196658_Sudhindra/comparative_genomics/ani_dir

average_nucleotide_identity.py -i ./ -o ./ANIm_output -m ANIm \
                               -g --gformat png,pdf
