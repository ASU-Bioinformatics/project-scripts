#!/bin/bash

##### QC of MAG bins with BUSCO #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.busco.out
#SBATCH -o slurm.%j.busco.err
#SBATCH -t 0-0:15

#### initiate array of busco checks ####
inputDir="/data/gencore/analysis_projects/8718394_Alsanea/metawrap-bins/Bucket/refined_bins/metawrap_70_5_bins"

list=$(find "$inputDir" -maxdepth 1 -type f -name "*.fa" | while read F; do basename $F; done | sort | uniq)
samples=($list)
count=${#samples[@]}

echo $list
echo $count

busco=$(sbatch --array=0-$((count-1)) --parsable --export=NONE \
        /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_G/G7.busco.sh "$list" "$inputDir")
