#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.B1.init.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.B1.init.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00
#SBATCH --mem=128G

module load r-4.2.2-gcc-11.2.0

help="FALSE"

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
    -h | --help)
        help="TRUE"
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

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script is used to run prewritten R scripts for differential gene expression analysis
  The current tools available are edgeR, DESeq2, and NOISeq.

  usage: bash B1.DEG.rscripts.sh -d /path/to/output/folder \
                                 -g /path/to/matrix.csv \
                                 -c /path/to/comparisons.csv \
                                 -s "deseq2 edger noiseq"

  options:
    [ -d  |   --directory   |   folder name for DE analysis output                                                 ]
    [ -g  |   --geneMatrix  |   full pathname (ideally) to gene count matrix csv file from stringtie               ]
    [ -c  |   --comparisons |   full pathname (ideally) to csv file containing the pairwise comparison groupings   ]
    [ -s  |   --script      |   a quoted list of all scripts to run, case insenstive                               ]
    [ -h  |   --help        |   prints informational message and exits                                             ]

EOF
exit;
fi

mkdir -p "$directory"

for i in $script
do
  echo $i
  if [[ $i == "deseq2" ]]; then
    Rscript /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B2.deseq2.R \
      -d $directory -g $geneMatrix -c $comparisons

  elif [[ $i == "edger" ]]; then
    Rscript /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B3.edgeR.R \
      -d $directory -g $geneMatrix -c $comparisons

  elif [[ $i == "noiseq" ]]; then
    Rscript /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B4.noiseq.R \
      -d $directory -g $geneMatrix -c $comparisons

  else
    echo "The script name provided cannot be processed."
  fi
done;
