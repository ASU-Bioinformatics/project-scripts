#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.G3.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.G3.err               # STDERR (%j = JobId)
#SBATCH -t 0-12:00                     # estimated time needed
#SBATCH --mem=64G

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/checkm-env

VALID_ARGS=$(getopt -o b:o: \
                    --long binDir:,outDir: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -b | --binDir)
        echo "Metabat bin fasta files will be read from '$2'"
        binDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "The output files will be written to '$2'"
        outDir="$2"
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


cd "$binDir"

checkm lineage_wf ./ ./checkm_lineage-test -x fa

mv ./checkm_lineage-test $outDir
