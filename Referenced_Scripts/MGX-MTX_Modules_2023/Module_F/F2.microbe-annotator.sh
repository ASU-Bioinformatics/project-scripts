#!/bin/bash

##### download and test microbe-annotator #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH -t 2-0:00
#SBATCH -c 16
#SBATCH --mem=120G

num=0"${SLURM_ARRAY_TASK_ID}"
echo "$num"

echo "microbe annotator for contig coassembly"
echo "trying lowmem settings: 6 days, 120G, 10 cores, 16 threads"
#echo "trying highmem settings: 2 days, 512G, 10 cores, 64 threads"

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/microbe-annotator-env/

#### now for annotation! ####

dbdir="$1"
proteinFasta="$2"
continue="$3"
outdir="$(dirname $proteinFasta)"/microbe-annotator-out-"$num"

mkdir -p "$outdir"

if [[ "$continue" == "TRUE" ]]
then
  microbeannotator -i ${proteinFasta%.faa}."$num".faa \
                   -d "$dbdir" -o "$outdir" -m diamond -p 1 -t 100 --continue_run

elif [[ "$continue" == "FALSE" ]]
then
  microbeannotator -i ${proteinFasta%.faa}."$num".faa \
                   -d "$dbdir" -o "$outdir" -m diamond -p 1 -t 100

fi
