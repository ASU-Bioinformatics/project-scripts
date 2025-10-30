#!/bin/bash

#SBATCH -o slurm.%j.B4.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.B4.err               # STDERR (%j = JobId)

### for sol ###
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH --mem=16G

krakenBin="/data/biocore/programs/KrakenTools-1.2"
kronaBin="/data/biocore/programs/KronaTools-2.8.1/bin"

VALID_ARGS=$(getopt -o i:o:h:k:n: \
                    --long inputDir:,outDir:,htmlDir:,krakenBin:,kronaBin: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Will look in the directory '$2' for braken report files"
        inputDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Tabular Krona results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -h | --htmlDir)
        echo "HTML Krona plots will be output to the directory '$2'"
        htmlDir="$2"
        shift 2
        ;;
    -k | --krakenBin)
        echo "Will use an installation of KrakenTools in the directory '$2'"
        krakenBin="$2"
        shift 2
        ;;
    -n | --kronaBin)
        echo "Will use an installation of KronaTools in the directory '$2'"
        kronaBin="$2"
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
  This script uses the Bracken species-level report files to visualize taxonomic population composition in interactive HTML files.

  usage: sbatch B4.krona-plots.sh
            --inputDir /path/to/bracken-reports
            --outDir /path/to/tabular-output
            --krakenBin /path/to/krakenTools/installation
            --kronaBin /path/to/kronaTools/installation
            (--help)

  options:
    [ -i  |   --inputDir  |   directory containing the Bracken report files (file name ending in S.breport)  ]
    [ -o  |   --outDir    |   output directory for Krona tabular files                                            ]
    [ -h  |   --htmlDir   |   output directory for Krona HTML files ]
    [ -k  |   --krakenBin |   directory for KrakenTools installation; default is /data/biocore/programs/KrakenTools-1.2       ]
    [ -n  |   --kronaBin  |   directory for KronaTools installation; default is /data/biocore/programs/KronaTools-2.8.1/bin ]
    [ -h  |   --help      |   prints an informational message and exits script                                                         ]
EOF
  exit;
fi

module load mamba/latest

mkdir -p "$outDir"
mkdir -p "$htmlDir"

for i in $(find "$inputDir" -maxdepth 1 -type f -name "*S.breport" | while read F; do basename $F; done | rev | cut -d "_" -f 2-10 | rev | sort | uniq)
do
  python "$krakenBin"/kreport2krona.py \
              -r "$inputDir"/"$i"_S.breport \
              -o "$outDir"/"$i"_S.b.krona.txt \
              --no-intermediate-ranks
  "$kronaBin"/ktImportText "$outDir"/"$i"_S.b.krona.txt \
              -o "$htmlDir"/"$i"_S.b.krona.html
done
