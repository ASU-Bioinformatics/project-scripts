#!/bin/bash

#SBATCH -o slurm.%j.B1.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.B1.err               # STDERR (%j = JobId)

### for sol ###
#SBATCH -p public
#SBATCH -q public
#SBATCH -t 1-00:00:00                 # took <1hr for a single sample, >50M usable reads and plusPF database
#SBATCH --mem=120G

# define project-specific variables

VALID_ARGS=$(getopt -o i:o:d: \
                    --long inputDir:,outDir:,dbDir: \
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
    -d | --dbDir)
        echo "The KRAKEN database in the directory '$2' will be used for classification"
        dbDir="$2"
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
  This script removes adapter sequences and filters for quality in 18S, 16S, and ITS reads. Honestly, it could be used
  for a variety of read types beyond that. The cutadapt parameters can't be customized at this time, and are -m 30 -q 20

  usage: sbatch B1.kraken.sh
            --inputDir /path/to/fastq-input
            -outDir /path/to/kraken-output
            -dbDir /path/to/kraken-db (-h)

  options:
    [ -i  |   --inputDir  |   directory containing fastq.gz files with traditional Illumina/Casava naming conventions                  ]
    [ -o  |   --outDir    |   directory for kraken output folders (raw output, reports, classified fastq, and unclassified fastq)      ]
    [ -d  |   --dbDir     |   directory containing the Kraken database of choice; default is /data/gencore/databases/kraken/k2_pluspf/ ]
    [ -h  |   --help      |   prints an informational message and exits script                                                         ]
EOF
  exit;
fi


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
  /data/biocore/programs/kraken2/kraken2 --db "$dbDir" --threads 8 \
                  --report "$reportDir"/"$sid"_SQP_k2pf_pe.k2report \
                  --report-minimizer-data --minimum-hit-groups 3 --paired \
                  --classified-out "$clsfDir"/"$sid"_k2pf_pe_clsf#.fq \
                  --unclassified-out "$unclDir"/"$sid"_SQP_k2pf_pe_uncl#.fq \
                  "$inputDir"/"$sid"_SQP_L00*_R1_001.fastq.gz \
                  "$inputDir"/"$sid"_SQP_L00*_R2_001.fastq.gz > "$rawDir"/"$sid"_SQP_k2pf_pe.kraken2
done;
