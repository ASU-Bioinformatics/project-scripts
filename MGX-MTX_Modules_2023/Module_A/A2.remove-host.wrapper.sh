#!/bin/bash

##### alignment to host reference #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-0:15

#### notes on A2.bowtie-align ####
# this script takes a long time to run, will need to be run for each sample independently
# set $assembly to the directory containing the bowtie2 index, including the index prefix
# set $fqDir to the directory containing the paired end fastq files, preferably quality trimmed with adapters removed
# set $outDir to the output directory for the sample alignments
# $use lets the wrapper script know whether to include all the samples in $fqDir (DIRECTORY) or not (LIST)
# $list provides the list of variables to use if the LIST option is specified for $use

#### Define Run Variables ####

list=""

VALID_ARGS=$(getopt -o a:f:o:u:l: \
                    --long assembly:,fqDir:,outDir:,use:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --assembly)
        echo "Assembly location and prefix will be found at '$2'"
        assembly="$2"
        shift 2
        ;;
    -f | --fqDir)
        echo "Input samples will be taken from the directory '$2'"
        fqDir="$2"
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

#assembly="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/megahit-alignments/coassembly"
#fqDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/fastq/adapter-trimmed/fastq"
#outDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/megahit-alignments/samples"
#use="DIRECTORY"
#list="56 58 62"

mkdir -p "$outDir"
mkdir -p "$outDir"/qc
mkdir -p "$outDir"/alignments

#### Call Script ####

if [ "$use" == "DIRECTORY" ];
then

  ### run A2.remove-host.sh for all samples in a directory ###

  for i in $(find "$fqDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".A2.remove-host \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A2.remove-host.sh \
         "$i" "$assembly" "$fqDir" "$outDir"
  done;

elif [ "$use" == "LIST" ];
then

  ### run A2.remove-host.sh for a subset of samples in a directory ###

  for i in $list
  do
    echo "$i"
    sbatch --job-name "$i".A2.remove-host \
         /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A2.remove-host.sh \
         "$i" "$assembly" "$fqDir" "$outDir"
  done;

fi
