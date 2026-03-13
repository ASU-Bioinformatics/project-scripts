#!/bin/bash

##### classify MAG bins with BAT #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH -o slurm.%j.G4.out
#SBATCH -o slurm.%j.G4.err
#SBATCH -t 4-0:00
#SBATCH -c 16
#SBATCH --mem=1024G

#### run CAT contigs on final coassembly fasta file ####

# updates from 2023 version include installing a new version of CAT as well as creating a new nr database as of 2023-03-11
# I also added arguments to specify the environment and script directory
# and created a help method to display on the command line

module load mamba/latest

catEnv="/data/biocore/programs/mamba-envs/catbat-env"
scriptDir="/data/biocore/programs/CAT_pack-6.0.1/CAT_pack"

# binDir is the absolute path to the directory containing the bin files with extension .fa
# databaseDir is the absolute path to the CAT database directory
# taxDir is the absolute path to the CAT taxonomy directory
# outDir is the absolute path to the directory where BAT output files should end up (defaults to current directory)

##### Define Variables #####

VALID_ARGS=$(getopt -o b:d:t:o:s:e:h \
                    --long binDir:,databaseDir:,taxDir:,outDir:,scriptDir:,catEnv:,help \
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
    -s | --scriptDir)
        echo "The CAT_pack scripts are in the folder '$2'"
        scriptDir="$2"
        shift 2
        ;;
    -e | --catEnv)
        echo "The mamba environment containing the CAT_pack dependences is '$2'"
        mwEnv="$2"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
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

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script uses metaWRAP to create metagenomic bins from an assembly fasta file using one of three different binning tools.
  To save time, only one tool can be called at a time; the options are 'concoct', 'metabat', and 'maxbin'.

  CheckM2 is run on the bins afterwards to predict the completeness and contamination of each one.

  usage: sbatch G3.bat-bin-annotate.sh
            -b /path/to/bin/directory
            -d /path/to/database/directory
            -t /path/to/taxonomy/directory
            -o /path/to/output/directory
            (-s /path/to/script/directory) (-e /path/to/catpack/environment) (-h)

  options:
    [ -b  | --binDir      | pathname for the directory containing fasta bins (suffix '.fa')                                               ]
    [ -d  | --databaseDir | pathname for the directory containing the prepared CAT database files                                         ]
    [ -t  | --taxDir      | pathname for the directory containing the prepared CAT taxonomy files                                         ]
    [ -o  | --outDir      | directory for output files (classification and protein predictions)                                           ]
    [ -s  | --scriptDir   | directory containing the CAT_pack scripts; default "/data/biocore/programs/CAT_pack-6.0.1/CAT_pack/CAT_pack"  ]
    [ -e  | --catEnv      | metaWRAP environment; default "/data/biocore/programs/mamba-envs/metawrap-env"                                ]
    [ -h  | --help        | prints an informational message and exits script                                                              ]
EOF
  exit;
fi

source activate $catEnv

mkdir -p "$outDir"
cd "$outDir"

"$scriptDir"/CAT_pack bins \
  -b "$binDir" \
  -d "$databaseDir" \
  -t "$taxDir" \
  -s .fa --verbose --index_chunks 1

# add names to classification file and summarize

"$scriptDir"/CAT_pack add_names \
  -i "$outDir"/out.BAT.bin2classification.txt \
  -o "$outDir"/out.BAT.bin2classification.official.txt \
  -t "$taxDir" --only_official

"$scriptDir"/CAT_pack summarise \
  -i "$outDir"/out.BAT.bin2classification.official.txt \
  -o "$outDir"/out.BAT.summarize.txt
