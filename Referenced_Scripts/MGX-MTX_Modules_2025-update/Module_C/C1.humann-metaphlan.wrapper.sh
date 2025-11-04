#!/bin/bash

##### humann3 test #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-12:00

#### notes on C1.humann-metaphlan ####
# this script takes a long time to run, will need to be run for each sample independently
# set $concatenate to TRUE if the fastq files need to be concatenated still
# set $concatenate to FALSE if the fastq files have already been concatenated
# set $type to DNA for metagenomics or metatranscriptomics without paired metagenomics data
# set $type to RNA for metatranscriptomics with a paired metagenomics humann taxonomic profile
# set $database to u50 for uniref50 analysis
# set $database to u90 for uniref90 analysis
# !! this will set the database for the whole system, not just this run !!
# !! don't try to run two sets of analysis on different databases simultaneously !!
# use search mode uniref90 for well-characterized microbiomes (ie, human gut)
# if the total reads unmapped is low or the microbiome is uncharacterized, try uniref50
# $refDir contains the paired metagenomics taxonomic profiles for metatranscriptomics analysis; can provide "" for DNA runs
# $inputDir contains the un-concatenated paired end fastq files, preferably quality trimmed with adapters removed
# $use lets the wrapper script know whether to include all the samples in $inputDir (DIRECTORY) or not (LIST)
# $list provides the list of variables to use if the LIST option is specified for $use

#### Define Run Variables
# defaults
concatenate="TRUE"
type="DNA"
refDir=""
use="DIRECTORY"
resume="FALSE"

VALID_ARGS=$(getopt -o i:cpo:r:d:ml: \
                    --long inputDir:,concatenated,paired,outDir:,dnaRef:,dbType:,resume,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Fastq files will be taken from the directory '$2'"
        inputDir="$2"
        shift 2
        ;;
    -c | --concatenated)
        echo "Fastq files are already concatenated"
        concatenate="FALSE"
        shift
        ;;
    -p | --paired)
        echo "Data is from a paired metatranscriptomics sample set"
        type="RNA"
        shift
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -r | --dnaRef)
        echo "Paired DNA metaphlan output is in the directory '$2'"
        refDir="$2"
        shift 2
        ;;
    -d | --dbType)
        echo "The HUMAnN database with UniRef set '$2' will be used for classification"
        database="$2"
        shift 2
        ;;
    -m | --resume)
        echo "Will begin Humann/Metaphlan from a previous unfinished run"
        resume="TRUE"
        shift
        ;;
    -l | --list)
        echo "Samples from the list '$2' will be classified"
        list="$2"
        use="LIST"
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

module load mamba/latest
source activate /data/biocore/programs/conda-envs/humann-env/

if [[ "$database" == "u50" ]]; then
  humann_config --update database_folders protein /data/gencore/databases/humann/uniref/uniref
  mode="uniref50"
elif [[ "$database" == "u90" ]]; then
  humann_config --update database_folders protein /data/gencore/databases/humann/uniref90/
  mode="uniref90"
fi

if [ "$use" == "DIRECTORY" ];
then

  ### run C1.humann-metaphlan for all samples in a directory ###

  for sid in $(find "$inputDir" -type f -name "*.fastq")
  do
    gzip "$sid"
  done

  for i in $(find "$inputDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".C1.humann-metaphlan \
         /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_C/C1.humann-metaphlan.sh \
         "$i" "$refDir" "$outDir" "$inputDir" "$concatenate" "$type" "$mode" "$resume"
  done;

elif [ "$use" == "LIST" ];
then

  ### to run C1.humann-metaphlan for a subset of samples in directory ###

  for i in $list
  do
    echo "$i"
    sbatch --job-name "$i".C1.humann-metaphlan \
         /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_C/C1.humann-metaphlan.sh \
         "$i" "$refDir" "$outDir" "$inputDir" "$concatenate" "$type" "$mode" "$resume"
  done;

fi
