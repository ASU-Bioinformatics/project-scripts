#!/bin/bash

##### QC of MAG bins with BUSCOclassify MAG bins with CheckM2 #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.busco.out
#SBATCH -o slurm.%j.busco.err
#SBATCH -t 0-2:00

#### check quality of bins using busco ####

echo $SLURM_ARRAY_TASK_ID
echo $1
samples=($1)
sid=${samples[$SLURM_ARRAY_TASK_ID]}

echo $samples

echo $sid

inDir=$2
cd $inDir

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/busco-env

busco -i $sid --auto-lineage -o ${sid%.fa} -m geno
