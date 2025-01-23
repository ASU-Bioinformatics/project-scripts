#!/bin/bash

##### alignment to host reference #####

###SBATCH -p htc
###SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-0:15

#### notes on A3.remove-host-rna ####
# this script takes a long time to run, will need to be run for each sample independently
# set $refDir to the directory containing the STAR indexes
# set $inputDir to the directory containing the paired end fastq files, quality trimmed with adapters removed
# set $outDir to the output directory for the sample alignments
# set $tempDir to a data drive that can handle STAR I/O issues (on Sol, use scratch)
# $use lets the wrapper script know whether to include all the samples in $inputDir (DIRECTORY) or not (LIST)
# $list provides the list of sample IDs to use if the LIST option is specified for $use

# the default STAR version on Sol is currently 2.7.10a so I don't load any module or environment
# just make sure to use 2.7.10a indexes

#### Define Run Variables ####

list=""

VALID_ARGS=$(getopt -o i:r:t:o:u:l: \
                    --long inputDir:,refDir:,tempDir:,outDir:,use:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Fastq files will be obtained from the directory '$2'"
        inputDir="$2"
        shift 2
        ;;
    -r | --refDir)
        echo "The STAR indexes will be obtained from the directory '$2'"
        refDir="$2"
        shift 2
        ;;
    -t | --tempDir)
        echo "STAR analysis will be performed in the directory '$2'"
        tempDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -u | --use)
        echo "Fastq files will come from the provided '$2'"
        use="$2"
        shift 2
        ;;
    -l | --list)
        echo "Fastq sample IDs from '$2' will be analyzed"
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
mkdir -p "$outDir"/qc
mkdir -p "$outDir"/alignments

#### Call Script ####

if [ "$use" == "DIRECTORY" ];
then

  ### run A3.remove-host-rna.sh for all samples in a directory ###

  for i in $(find "$inputDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".A3.remove-host-rna \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A3.remove-host-rna.sh \
         "$i" "$inputDir" "$refDir" "$tempDir" "$outDir"
  done;

elif [ "$use" == "LIST" ];
then

  ### run A2.remove-host.sh for a subset of samples in a directory ###

  for i in $list
  do
    echo "$i"
    sbatch --job-name "$i".A3.remove-host-rna \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A3.remove-host-rna.sh \
         "$i" "$inputDir" "$refDir" "$tempDir" "$outDir"
  done;

fi
