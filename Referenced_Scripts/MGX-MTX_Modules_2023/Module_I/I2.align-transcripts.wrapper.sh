#!/bin/bash

##### align metatranscriptomic reads to paired metagenomic predicted genes #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err               # STDERR (%j = JobId)
#SBATCH -t 0-00:15

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/nanopore-env

# transcriptDir contains the filtered and trimmed data
# outDir is the folder where all output files are stored
# refDir is the folder containing the indexes
# ref is the index ID, along with the absolute pathname from refDir
# use is either DIRECTORY or LIST, depending on whether all files or only some files from a directory need to be aligned
# list provides the list of samples if LIST is chosen

# default values
use="DIRECTORY"
list=""

VALID_ARGS=$(getopt -o t:o:h:l: \
                    --long transcriptDir:,outDir:,hisatIndex:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -t | --transcriptDir)
        echo "Fastq reads from paired transcriptomes will be taken from the directory '$2'"
        transcriptDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -h | --hisatIndex)
        echo "The hisat index with the pathname and prefix '$2' will be used."
        hisatIndex="$2"
        shift 2
        ;;
    -l | --list)
        echo "Samples will be taken from the following list of IDs present in the transcript directory: '$2'"
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

mkdir -p "$outDir"

if [[ $use == "DIRECTORY" ]]
then

  #### run I2.align-transcripts.sh for all samples in $transcriptDir ####

  for i in $(find "$transcriptDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".I2.align-transcripts \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I2.align-transcripts.sh \
         "$transcriptDir" "$hisatIndex" "$i" "$outDir"
  done;

elif [[ $use == "LIST" ]]
then

  #### run I2.align-transcripts.sh only on samples specified in $list

  for i in "$list"
  do
    echo "$i"
    sbatch --job-name "$i".I2.align-transcripts \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I2.align-transcripts.sh \
         "$transcriptDir" "$hisatIndex" "$i" "$outDir"
  done;

else

  echo "variable 'use' must be either LIST or DIRECTORY"

fi
