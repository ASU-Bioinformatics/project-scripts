#!/bin/bash

##### refine bins from three binners run with G1 #####

#SBATCH -p public
#SBATCH -q public
#SBATCH --output slurm.%j.MGX.G2.out
#SBATCH --error slurm.%j.MGX.G2.err
#SBATCH -c 24
#SBATCH -t 2-0:00
#SBATCH --mem=120G
###SBATCH -q debug
###SBATCH -p htc
###SBATCH -t 0-0:15

umask 0007
module purge
module load mamba/latest

##### Define Variables #####

mwEnv="/data/biocore/programs/mamba-envs/metawrap-env"
minComp=70
maxCont=5

VALID_ARGS=$(getopt -o a:b:c:o:m:x:e:h \
                    --long binDirA:,binDirB:,binDirC:,outDir:,minComp:,maxCont:,metawrapEnvironment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --binDirA)
        echo "The first set of bins to include in refinement is found in '$2'"
        binDirA="$2"
        shift 2
        ;;
    -b | --binDirB)
        echo "The second set of bins to include in refinement is found in '$2'"
        binDirB="$2"
        shift 2
        ;;
    -c | --binDirC)
        echo "The third set of bins to include in refinement is found in '$2'"
        binDirC="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Bin refinement output will be written to '$2'"
        outDir="$2"
        shift 2
        ;;
    -m | --minComp)
        echo "The minimum completion percentage is '$2'"
        minComp="$2"
        shift 2
        ;;
    -x | --maxCont)
        echo "The maximum contamination percentage is '$2'"
        maxCont="$2"
        shift 2
        ;;
    -e | --metawrapEnvironment)
        echo "The conda environment to use for metaWRAP is '$2'"
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

  usage: sbatch G2.metaWRAP-refinement.sh
            -a /path/to/binnerA/bins
            -b /path/to/binnerB/bins
            -c /path/to/binnerC/bins
            -o /path/to/output/directory
            -m minimumCompletionPercentage -x maximumContaminationPercentage
            (-e /path/to/metawrap/environment) (-h)

  options:
    [ -a  | --binDirA             | pathname for the directory containing identified bins from tool A               ]
    [ -b  | --binDirB             | pathname for the directory containing identified bins from tool B               ]
    [ -c  | --binDirC             | pathname for the directory containing identified bins from tool C               ]
    [ -o  | --outDir              | directory for output files (stats, refined bins, etc.)                          ]
    [ -m  | --minComp             | minimum completion percentage for consolidated bin (default is 70%)             ]
    [ -x  | --maxCont             | maximum contamination percentage for consolidated bin (default is 5%)           ]
    [ -e  | --metawrapEnvironment | metaWRAP environment; default "/data/biocore/programs/mamba-envs/metawrap-env"  ]
    [ -h  | --help                | prints an informational message and exits script                                ]
EOF
  exit;
fi

#### run metaWrap binning tool with specified binner #####

source activate "$mwEnv"

metawrap bin_refinement \
  -o $outDir -t 96 -m 360 \
  -c $minComp -x $maxCont \
  -A $binDirA -B $binDirB -C $binDirC

source deactivate
