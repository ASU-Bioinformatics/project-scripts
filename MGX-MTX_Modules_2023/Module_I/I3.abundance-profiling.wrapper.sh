#!/bin/bash

##### I3.abundance-profiling.sh wrapper script #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-0:15

#### notes on I3.abundance-profiling ####
# this script takes a long time to run, will need to be run for each sample independently (hence the wrapper)
# set $alnDir to the directory containing the SAM files from each sample's alignment to the coassembly fasta file
# set $gffFile to the CAT predicted protein gff file corresponding to the coassembly fasta file
# set $outDir to the directory into which the final htseq count files should be placed
# $use lets the wrapper script know whether to include all the samples in $alnDir (DIRECTORY) or not (LIST)
# $list provides the list of variables to use if the LIST option is specified for $use

use="DIRECTORY"
list=""

VALID_ARGS=$(getopt -o a:o:g:l: \
                    --long alignmentDir:,outDir:,gffFile:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --alignmentDir)
        echo "Transcript alignments will be read from '$2'"
        alnDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -g | --gffFile)
        echo "Coassembly annotations will be read from '$2'"
        gffFile="$2"
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

  #### run I3.abundance-profiling.sh for all samples in $alnDir ####

  for i in $(find "$alnDir" -type f -name "*.sorted.bam" | while read F; do basename $F; done | rev | cut -d "." -f 3-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".I3.abundance-profiling \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I3.abundance-profiling.sh \
         "$i" "$alnDir" "$gffFile" "$outDir"
  done;

elif [[ $use == "LIST" ]]
then

  #### run I3.abundance-profiling.sh only on samples specified in $list

  for i in "$list"
  do
    echo "$i"
    sbatch --job-name "$i".I3.abundance-profiling \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I3.abundance-profiling.sh \
         "$i" "$alnDir" "$gffFile" "$outDir"
  done;

else

  echo "variable 'use' must be either LIST or DIRECTORY"

fi
