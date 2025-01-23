#!/bin/bash

##### download and test interproscan #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH -o slurm.%j.G2.out
#SBATCH -o slurm.%j.G2.err
#SBATCH -t 4-0:00
#SBATCH -c 16
#SBATCH --mem=1024G

#### run CAT contigs on final coassembly fasta file ####

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/catbat-env

# binDir is the absolute path to the directory containing the bin files with extension .fa
# databaseDir is the absolute path to the CAT database directory
# taxDir is the absolute path to the CAT taxonomy directory
# outDir is the absolute path to the directory where BAT output files should end up (defaults to current directory)

##### Define Variables #####

VALID_ARGS=$(getopt -o b:d:t:o: \
                    --long binDir:,databaseDir:,taxDir:,outDir: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -b | --binDir)
        echo "Metabat bin fasta files will be read from '$2'"
        binDir="$2"
        shift 2
        ;;
    -d | --databaseDir)
        echo "The CAT database for protein prediction is '$2'"
        databaseDir="$2"
        shift 2
        ;;
    -t | --taxDir)
        echo "The CAT taxonomy for contig classification is '$2'"
        taxDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "The output files will be written to '$2'"
        outDir="$2"
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
cd "$outDir"

CAT bins -b "$binDir" \
         -d "$databaseDir" \
         -t "$taxDir" \
         -s .fa --verbose --index_chunks 1

# add names to classification file and summarize

CAT add_names -i "$outDir"/out.BAT.bin2classification.txt \
              -o "$outDir"/out.BAT.bin2classification.official.txt \
              -t "$taxDir" --only_official

CAT summarise -i "$outDir"/out.BAT.bin2classification.official.txt \
              -o "$outDir"/out.BAT.summarize.txt
