#!/bin/bash

##### download and test microbe-annotator #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.F2-split.out
#SBATCH -e slurm.%j.F2-split.err
#SBATCH -t 0-00:15
#SBATCH -c 1
#SBATCH --mem=16G

# set $dbdir to the directory containing the microbe-annotator database
# set $fadir to the directory containing the predicted proteins .faa file from the cat contig annotator
# $sDir should be the directory containing the F3.microbe-annotator script
# set $continue to TRUE if the run is being resumed after a crash or timeout error
# set $continue to FALSE if beginning a run for the first time
# set $split to TRUE if you need to split the predicted proteins file into smaller chunks
# set $split to FALSE if the file has already been split and you're running or rerunning portions of it
# if $split is FALSE, indicate the chunk indexes with the variable $list, which should be comma separated

echo "microbe annotator for contig coassembly - split file into smaller chunks if needed"

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/microbe-annotator-env/

#### Define Variables

dbdir="/data/gencore/databases/microbe-annotator"
sDir="/data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F"
split="FALSE"
continue="FALSE"

use="DIRECTORY"
list=""

VALID_ARGS=$(getopt -o d:p:scl: \
                    --long databaseDir:,proteinFasta:,split,continue,list \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -d | --databaseDir)
        echo "MicrobeAnnotator database is located at '$2'."
        dbDir="$2"
        shift 2
        ;;
    -p | --proteinFasta)
        echo "Predicted proteins to be annotated are found at '$2'."
        proteinFasta="$2"
        shift 2
        ;;
    -s | --split)
        echo "Coassembly file will be split into 10 sections for faster processing."
        split="TRUE"
        shift
        ;;
    -c | --continue)
        echo "Annotation will continue an incomplete run."
        continue="TRUE"
        shift
        ;;
    -l | --list)
        echo "Samples will be taken from the following list of split assembly sections present in the transcript directory: '$2'"
        use=LIST
        list="$2"
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


if [[ "$split" == "TRUE" ]]
then

  #### split the coassembly if it is too large ####

  cd "$fadir"
  pyfasta split -n 10 "$proteinFasta"

  #### now for annotation! ####

  sbatch --array=0-9%5 --parsable \
       -p general -q public --job-name F2.microbe-annotator \
       --export=ALL "$sDir"/F2.microbe-annotator.sh \
       "$dbdir" "$proteinFasta" "$continue"

elif [[ "$split" == "FALSE" ]]
then

  if [[ "$use" == "LIST" ]]
  then

    #### use the list to signify the array indexes to run ####

    sbatch --array="$list" --parsable \
       -p general -q public --job-name F2.microbe-annotator \
       --export=ALL "$sDir"/F2.microbe-annotator.sh \
       "$dbdir" "$proteinFasta" "$continue"

  elif [[ "$use" == "DIRECTORY" ]]
  then

    #### use all files 0-9 from the split ####

    sbatch --array=0-9%5 --parsable \
         -p general -q public --job-name F2.microbe-annotator \
         --export=ALL "$sDir"/F2.microbe-annotator.sh \
         "$dbdir" "$proteinFasta" "$continue"
  fi
fi
