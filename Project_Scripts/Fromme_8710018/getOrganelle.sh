#!/bin/bash

##### getOrganelle assembler #####

#SBATCH -p general
####SBATCH -q grp_kawoodbu
#SBATCH -q public
###SBATCH -q debug
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
###SBATCH -t 0-0:10
#SBATCH -t 0-12:00                         # estimated time needed - edited down from 2-0:00 after first run
#SBATCH --mem=32G                        # again potentially overloading this - edited down from 256G after first run
#SBATCH -c 1                              # I don't think I'll need two CPUs

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/getOrganelle-env

# upload databases
# python /data/biocore/programs/mamba-envs/getOrganelle-env/bin/get_organelle_config.py -a all --config-dir /data/gencore/databases/getOrganelle/getOrganelle

# test
# get_organelle_from_reads.py -1 Arabidopsis_simulated.1.fq.gz -2 Arabidopsis_simulated.2.fq.gz \
#  -t 1 -o Arabidopsis_simulated.plastome -F embplant_pt -R 10 --config-dir /data/gencore/databases/getOrganelle/getOrganelle

#get_organelle_from_reads.py \
#  -1 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R1_001.fastq.gz \
#  -2 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R2_001.fastq.gz \
#  -o /data/gencore/analysis_projects/8710018_Sunidhi/getOrganelle_output \
#  -R 15 -k 21,45,65,85,105 -F embplant_pt --config-dir /data/gencore/databases/getOrganelle/getOrganelle

# changes: no set max rounds (-R); allowing default based on (-F) param
# include embplant_mt seed database
#get_organelle_from_reads.py \
#  -1 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R1_001.fastq.gz \
#  -2 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R2_001.fastq.gz \
#  -o /data/gencore/analysis_projects/8710018_Sunidhi/getOrganelle_output_params2 \
#  -k 17,37,57,77,97,117 -F embplant_pt,embplant_mt --config-dir /data/gencore/databases/getOrganelle/getOrganelle

# changes - same as round 2 but with different kmer lengths
get_organelle_from_reads.py \
  -1 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R1_001.fastq.gz \
  -2 /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/Chloroplast_SCT_L001_R2_001.fastq.gz \
  -o /data/gencore/analysis_projects/8710018_Sunidhi/getOrganelle_output_params3 \
  -k 19,35,51,69,83 -F embplant_pt,embplant_mt --config-dir /data/gencore/databases/getOrganelle/getOrganelle
