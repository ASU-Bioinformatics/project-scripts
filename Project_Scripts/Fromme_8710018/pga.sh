#!/bin/bash

#### fastqc file generation #####

#SBATCH -p general
#SBATCH -q grp_kawoodbu #phx
#SBATCH -t 2-0:00               # estimated time needed
#SBATCH --mem=256G
#SBATCH -c 2

umask 0007

module load blast-plus-2.12.0-6z

perl /data/biocore/programs/PGA/PGA.pl \
	-r /data/gencore/analysis_projects/8710018_Sunidhi/cactaceae-chloroplasts \
	-t /data/gencore/analysis_projects/8710018_Sunidhi/getOrganelle_output/pga
