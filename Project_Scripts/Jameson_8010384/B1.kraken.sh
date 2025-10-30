#!/bin/bash

#SBATCH -o slurm.%j.B1.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.B1.err               # STDERR (%j = JobId)

### for sol ###
#SBATCH -p general
#SBATCH -q public
#SBATCH -t 1-00:00:00                 # took <1hr for a single sample, >50M usable reads and plusPF database
#SBATCH --mem=120G

# define project-specific variables

VALID_ARGS=$(getopt -o i:o:d:h: \
                    --long inputDir:,outDir:,dbDir:,help: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

help="FALSE"

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
    -d | --dbDir)
        echo "The KRAKEN database in the directory '$2' will be used for classification"
        dbDir="$2"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
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

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script is used to run Kraken assembly-free metagenomic analysis on all fastq files in the specified directory.

  usage: sbatch ./B1.kraken.sh \
          --inputDir /pathway/to/fastq \
          --outDir /pathway/to/kraken-output \
          --dbDir /data/gencore/databases/kraken/k2_pluspf

  options:
    [ -i  |   --inputDir  |   pathway to folder containing input fastq.gz files                             ]
    [ -o  |   --outDir    |   pathway to folder in which Kraken output files should be written              ]
    [ -d  |   --dbDir     |   pathway to the folder containing the desired Kraken classification database   ]
    [ -h  |   --help      |   prints informational message and exits                                        ]

EOF
exit;
fi

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/kraken2-env-2025

rawDir="$outDir"/kraken-raws
reportDir="$outDir"/kraken-reports
clsfDir="$outDir"/classified-fastq
unclDir="$outDir"/unclassified-fastq

mkdir -p $rawDir
mkdir -p $reportDir
mkdir -p $clsfDir
mkdir -p $unclDir

cd "$inputDir"


for sid in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
do
  echo "$sid"
  kraken2 --db "$dbDir" --threads 8 \
          --report "$reportDir"/"$sid"_SQP_k2pf_pe.k2report \
          --report-minimizer-data --minimum-hit-groups 3 --paired \
          --classified-out "$clsfDir"/"$sid"_k2pf_pe_clsf#.fq \
          --unclassified-out "$unclDir"/"$sid"_SQP_k2pf_pe_uncl#.fq \
          "$inputDir"/"$sid"_SQP_L00*_R1_001.fastq.gz \
          "$inputDir"/"$sid"_SQP_L00*_R2_001.fastq.gz > "$rawDir"/"$sid"_SQP_k2pf_pe.kraken2
done;
