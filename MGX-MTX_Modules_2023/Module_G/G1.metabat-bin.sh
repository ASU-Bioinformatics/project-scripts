#!/bin/bash

##### metabat binning #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.G1.out
#SBATCH -e slurm.%j.G1.err
#SBATCH -c 4
#SBATCH -t 3-0:00
#SBATCH --mem=128G

module purge
module load mamba/latest

source activate /data/biocore/programs/conda-envs/megahit-args

##### Define Variables #####

VALID_ARGS=$(getopt -o b:c: \
                    --long bamDir:,contigs: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -b | --bamDir)
        echo "Alignment sorted bam files will be read from '$2'"
        bamDir="$2"
        shift 2
        ;;
    -c | --contigs)
        echo "The contig assembly for binning is '$2'"
        contigs="$2"
        shift 2
        ;;
    --)
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1 - please correct."
        ;;
  esac
done

samples=$(find "$bamDir" -maxdepth 1 -type f -name "*.sorted.bam" | while read F; do basename $F; done)

cd "$bamDir"

runMetaBat.sh "$contigs" \
              $samples
