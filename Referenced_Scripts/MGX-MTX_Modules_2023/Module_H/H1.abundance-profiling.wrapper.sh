#!/bin/bash

##### H1.abundance-profiling.sh wrapper script #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-0:15

#### notes on H1.abundance-profiling ####
# this script takes a long time to run, will need to be run for each sample independently (hence the wrapper)
# set $bamDir to the directory containing the SAM files from each sample's alignment to the coassembly fasta file
# set $gffFile to the CAT predicted protein gff file corresponding to the coassembly fasta file
# set $outDir to the directory into which the final htseq count files should be placed
# $use lets the wrapper script know whether to include all the samples in $bamDir (DIRECTORY) or not (LIST)
# $list provides the list of variables to use if the LIST option is specified for $use

use="DIRECTORY"
list=""

VALID_ARGS=$(getopt -o b:g:o:l \
                    --long bamDir:,gffFile:,outDir:,list \
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
    -g | --gffFile)
        echo "The GFF file containing the predicted proteins is '$2'"
        gffFile="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Output files will be written to '$2'"
        outDir="$2"
        shift 2
        ;;
    -l | --list)
        echo "The list of files to be counted from the bam directory is '$2'"
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


mkdir -p "$outDir"

if [[ $use == "DIRECTORY" ]]
then

  #### run H1.abundance-profiling.sh for all samples in $bamDir ####

  for i in $(find "$bamDir" -type f -name "*.sorted.bam" | while read F; do basename $F; done | rev | cut -d "." -f 3-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".H1.abundance-profiling \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_H/H1.abundance-profiling.sh \
         "$i" "$bamDir" "$gffFile" "$outDir"
  done;

elif [[ $use == "LIST" ]]
then

  #### run H1.abundance-profiling.sh only on samples specified in $list

  for i in "$list"
  do
    echo "$i"
    sbatch --job-name "$i".H1.abundance-profiling \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_H/H1.abundance-profiling.sh \
         "$i" "$bamDir" "$gffFile" "$outDir"
  done;

else

  echo "variable 'use' must be either LIST or DIRECTORY"

fi
