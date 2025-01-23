#!/bin/bash

##### humann3 test #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-0:15

#### notes on D3.bowtie-align ####
# this script takes a long time to run, will need to be run for each sample independently
# set $assembly to the full path name through the prefix used to create the bowtie2 index for the assembly
# set $fastqDir to the directory containing the paired end fastq files, preferably quality trimmed with adapters removed
# set $outDir to the output directory for the sample alignments
# $list provides a list of variables to use (otherwise all samples in $fastqDir are used)

#### Define Run Variables ####

use="DIRECTORY"
list=""

VALID_ARGS=$(getopt -o f:a:o:l: \
                    --long fastqDir:,assembly:,outDir:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --fastqDir)
        echo "Fastq files will be taken from the directory '$2'"
        fqDir="$2"
        shift 2
        ;;
    -a | --assembly)
        echo "The assembly located at '$2'"
        assembly="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Alignments will be written to '$2'"
        outDir="$2"
        shift 2
        ;;
    -l | --list)
        echo "Samples will be taken from the list '$2' within the fastq directory"
        use="LIST"
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

#### Call Script ####

if [ "$use" == "DIRECTORY" ];
then

  ### run D3.bowtie-align for all samples in a directory ###

  for i in $(find "$fqDir" -maxdepth 1 -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".D3.bowtie-align \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D3.bowtie-align.sh \
         "$i" "$assembly" "$fqDir" "$outDir"
  done;

elif [ "$use" == "LIST" ];
then

  ### run D3.bowtie-align for a subset of samples in a directory ###

  for i in $list
  do
    echo "$i"
    sbatch --job-name "$i".D3.bowtie-align \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D3.bowtie-align.sh \
         "$i" "$assembly" "$fqDir" "$outDir"
  done;

fi
