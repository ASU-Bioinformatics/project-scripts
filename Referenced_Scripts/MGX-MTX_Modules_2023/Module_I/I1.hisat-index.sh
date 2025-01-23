#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.I1.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.I1.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00                       # for hisat only, can reduce to 2 days and 40G memory
#SBATCH -c 6
#SBATCH --mem=32G

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/biocore-rna

#### read in variables ####

VALID_ARGS=$(getopt -o o:p:c: \
                    --long outDir:,prefix:,contigs: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -o | --outDir)
        echo "Shortname FASTA file and HiSat2 index will be written to '$2'"
        outDir="$2"
        shift 2
        ;;
    -p | --prefix)
        echo "The prefix for the HiSat2 index will be '$2'"
        prefix="$2"
        shift 2
        ;;
    -c | --contigs)
        echo "The Hisat2 index will be made for the contig assembly file at '$2'"
        contigs="$2"
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

snFasta="$outDir"/shortname.final.contigs.fa
hsBase="$outDir"/"$prefix"

mkdir -p $outDir
cut -d " " -f1 $contigs > $snFasta

#### create index from predicted genes fasta file ####

### code to install HiSat2 if needed
# download and unpack HiSat2
#cd /data/biocore/programs
#wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download
#unzip download

# build index
# use --large-index tag if build failed with error requesting large index

/data/biocore/programs/hisat2-2.2.1/hisat2-build "$snFasta" "$hsBase" -p 64 --large-index
