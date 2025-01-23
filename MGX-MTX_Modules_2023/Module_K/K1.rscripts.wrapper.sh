#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.K1.rscripts.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.K1.rscripts..err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00
#SBATCH --mem=128G

module load r-4.2.2-gcc-11.2.0

## input variables from command line ##
VALID_ARGS=$(getopt -o d:g:c:s: \
                    --long directory:,geneMatrix:,comparisons:,script: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -d | --directory)
        echo "Output files will be written to '$2'"
        directory="$2"
        shift 2
        ;;
    -g | --geneMatrix)
        echo "The gene count matrix file can be found at '$2'"
        geneMatrix="$2"
        shift 2
        ;;
    -c | --comparisons)
        echo "The comparison csv file can be found at '$2'"
        comparisons="$2"
        shift 2
        ;;
    -s | --script)
        echo "The R script to call is '$2'"
        script="$2"
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

mkdir -p "$directory"

if [[ ${script,,} == "deseq2" ]]; then
  Rscript /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K2.deseq2.R \
          -d $directory -g $geneMatrix -c $comparisons

elif [[ ${script,,} == "edger" ]]; then
  Rscript /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K3.edgeR.R \
          -d $directory -g $geneMatrix -c $comparisons

elif [[ ${script,,} == "noiseq" ]]; then
  Rscript /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K4.noiseq.R \
          -d $directory -g $geneMatrix -c $comparisons

else
  echo "Please specify a script name from the following list: deseq2, edger, noiseq, clusterprofiler"
  break

fi
