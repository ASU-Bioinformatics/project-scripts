#!/bin/bash

##### ARG-focused metagenomics analysis using megahit #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH --output slurm.%j.D1.out              # STDOUT (%j = JobId)
#SBATCH -c 128
#SBATCH -t 2-0:00
#SBATCH --mem=1024G

module purge
module load mamba/latest

##### Define Variables #####

resume="FALSE"

VALID_ARGS=$(getopt -o i:o:r \
                    --long inputDir:,outDir:,resume \
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
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -r | --resume)
        echo "MEGAHIT will continue building the assembly in the output directory"
        resume="TRUE"
        shift
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

# project-specific variables


echo "running coassembly with MEGATHIT on samples in $inputDir, to output $outDir"

#### run megahit #####

source activate /data/biocore/programs/conda-envs/megahit-args

if [[ "$resume" = "" ]]; then
  forwards=$(find "$inputDir" -type f -name "*_SQP_L00*_R1_001.fastq.gz" | sort | uniq | tr '\n' ',')
  reverses=$(find "$inputDir" -type f -name "*_SQP_L00*_R2_001.fastq.gz" | sort | uniq | tr '\n' ',')
  megahit -1 $forwards \
          -2 $reverses \
          -o "$outDir" \
          --presets meta-large -t 128
else
  forwards=$(find "$inputDir" -type f -name "*_SQP_L00*_R1_001.fastq.gz" | sort | uniq | tr '\n' ',')
  reverses=$(find "$inputDir" -type f -name "*_SQP_L00*_R2_001.fastq.gz" | sort | uniq | tr '\n' ',')
  megahit -1 $forwards \
          -2 $reverses \
          -o "$outDir" --continue \
          --presets meta-large -t 128
fi

source deactivate
source activate /data/biocore/programs/conda-envs/quast-env

cd "$outDir"

quast -o coassembly-dna-quast final.contigs.fa
