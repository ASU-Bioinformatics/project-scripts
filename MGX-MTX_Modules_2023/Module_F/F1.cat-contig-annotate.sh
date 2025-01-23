#!/bin/bash

##### annotate assembled contigs with CAT #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH -o slurm.%j.F1.out
#SBATCH -o slurm.%j.F1.err
#SBATCH -t 0-8:00
#SBATCH -c 16
#SBATCH --mem=1024G

#### run CAT contigs on final coassembly fasta file ####

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/catbat-env

# default values
databaseDir="/data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database"
taxaDir="/data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy"

# commandline argument values
VALID_ARGS=$(getopt -o d:t:c:o: \
                    --long databaseDir:,taxaDir:,contigs:,outDir: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -d | --databaseDir)
        echo "CAT database is in the directory '$2'"
        databaseDir="$2"
        shift 2
        ;;
    -t | --taxaDir)
        echo "CAT taxa are in the directory '$2'"
        taxaDir="$2"
        shift 2
        ;;
    -c | --contigs)
        echo "The contig file '$2' will be used."
        contigs="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Output will be written to the directory '$2'"
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

#### project-specific variables ####

mkdir -p "$outDir"
cd "$outDir"

CAT contigs -c "$contigs" \
            -d "$databaseDir" \
            -t "$taxaDir" \
            --verbose --index_chunks 1 \
            --block_size 24 \
            --I_know_what_Im_doing --top 11 \
            -p out.CAT.predicted_proteins.faa \
            -a out.CAT.alignment.diamond

  # insert these tags if resuming from a checkpoint where they have been produced:

  #-p coassembly.CAT.predicted_proteins.faa \
  #-a coassembly.CAT.alignment.diamond \
