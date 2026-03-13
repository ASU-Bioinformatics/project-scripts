#!/bin/bash

##### refine bins from three binners run with G1 #####

#SBATCH -p public
#SBATCH -q public
#SBATCH --output slurm.%j.MGX.G3.reassemble.out
#SBATCH --error slurm.%j.MGX.G3.reassemble.err
#SBATCH -c 32
#SBATCH -t 2-0:00
#SBATCH --mem=196G
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

VALID_ARGS=$(getopt -o f:r:b:o:m:x:e:h \
                    --long forwards:,reverses:,bins:,outDir:,minComp:,maxCont:,metawrapEnvironment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --forwards)
        echo "The forward read fastq file is '$2'"
        forwards="$2"
        shift 2
        ;;
    -r | --reverses)
        echo "The reverse read fastq file is '$2'"
        reverses="$2"
        shift 2
        ;;
    -b | --bins)
        echo "The refined bins to reassemble are found in '$2'"
        bins="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Bin reassembly output will be written to '$2'"
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
  This script uses metaWRAP to reassemble the refined bins to improve assembly stats and reduce contamination.
  As with the metaWRAP binning step, the input reads need to be unzipped fastq with the naming format sample_1.fastq
  for the forward reads and sample_2.fastq for the reverse reads.

  CheckM2 is run on the bins afterwards to predict the completeness and contamination of each one.

  usage: sbatch G3.metaWRAP-reassemble.sh
            -f /path/to/fastq/sample_1.fastq
            -r /path/to/fastq/sample_2.fastq
            -b /path/to/refined/bins
            -o /path/to/output/directory
            -m minimumCompletionPercentage -x maximumContaminationPercentage
            (-e /path/to/metawrap/environment) (-h)

  options:
    [ -f  | --forwards            | pathname for the forward reads fastq file                                       ]
    [ -r  | --reverses            | pathname for the reverse reads fastq file                                       ]
    [ -b  | --bins                | pathname for the directory containing refined bins                              ]
    [ -o  | --outDir              | directory for output files (stats, reassembled bins, etc.)                      ]
    [ -m  | --minComp             | minimum completion percentage for consolidated bin (default is 90%)             ]
    [ -x  | --maxCont             | maximum contamination percentage for consolidated bin (default is 5%)           ]
    [ -e  | --metawrapEnvironment | metaWRAP environment; default "/data/biocore/programs/mamba-envs/metawrap-env"  ]
    [ -h  | --help                | prints an informational message and exits script                                ]
EOF
  exit;
fi

#### run metaWrap binning tool with specified binner #####

source activate "$mwEnv"

metawrap reassemble_bins \
  -o $outDir \
  -1 $forwards \
  -2 $reverses \
  -t 96 -m 180 \
  -c $minComp -x $maxCont \
  -b $bins

source deactivate
