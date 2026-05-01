#!/bin/bash

##### classify MAG bins with CheckM2 #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.checkm2.out
#SBATCH -o slurm.%j.checkm2.err
#SBATCH -t 4-0:00
#SBATCH -c 16
#SBATCH --mem=64G

#### classify and check quality of bins using checkm2 ####

# updates from 2023 version include installing a new version of CAT as well as creating a new nr database as of 2023-03-11
# I also added arguments to specify the environment and script directory
# and created a help method to display on the command line

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/checkm2-anaconda

checkm2 predict --threads 30 -x fa --force \
  --input /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/refined_bins/metawrap_70_5_bins \
  --output-directory /data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/refined_bins/checkm2_70_5_bins \
  --database_path /data/gencore/databases/CheckM2_database/uniref100.KO.1.dmnd
