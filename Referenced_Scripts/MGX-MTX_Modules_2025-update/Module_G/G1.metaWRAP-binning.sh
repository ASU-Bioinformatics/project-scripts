#!/bin/bash

##### metagenome assembly using megahit #####

#SBATCH -p public
#SBATCH -q public
#SBATCH --output slurm.%j.MGX.G1.out
#SBATCH --error slurm.%j.MGX.G1.err
#SBATCH -c 8
#SBATCH -t 1-0:00
#SBATCH --mem=196G
###SBATCH -q debug
###SBATCH -p htc
###SBATCH -t 0-0:15

module purge
module load mamba/latest

##### Define Variables #####

mwEnv="/data/biocore/programs/mamba-envs/metawrap-env"
chEnv="/data/biocore/programs/mamba-envs/checkm2-anaconda"
minLength=2000

VALID_ARGS=$(getopt -o f:r:o:a:b:e:c:l:h \
                    --long forwardReads:,reverseReads:,outputDirectory:,assembly:,binner:,metawrapEnvironment:,checkm2Environment:,minLength:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --forwardReads)
        echo "The forward reads are: '$2'"
        forwards="$2"
        shift 2
        ;;
    -r | --reverseReads)
        echo "The reverse reads are: '$2'"
        reverses="$2"
        shift 2
        ;;
    -o | --outputDirectory)
        echo "The bins will be output to '$2'"
        outDir="$2"
        shift 2
        ;;
    -a | --assembly)
        echo "The input assembly fasta file is '$2'"
        assembly="$2"
        shift 2
        ;;
    -b | --binner)
        echo "The binning tool to use is: '$2'"
        binner="$2"
        shift 2
        ;;
    -e | --metawrapEnvironment)
        echo "The conda environment to use for megaWRAP is '$2'"
        mwEnv="$2"
        shift 2
        ;;
    -c | --checkm2Environment)
        echo "The conda environment to use for checkM2 is '$2'"
        chEnv="$2"
        shift 2
        ;;
    -l | --minLength)
        echo "The minimum contig length to bin is: '$2'"
        minLength="$2"
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

  usage: sbatch G1.metaWRAP-binning.sh
            -f sample.r1.fastq.gz
            -r sample.r2.fastq.gz
            -o /path/to/output/directory
            -a assembly.fa
            -b "metabat" OR "maxbin" OR "concoct"
            -l 2000
            (-e /path/to/metawrap/environment) (-h)

  options:
    [ -f  | --forwardReads        | pathname for the original forward fastq reads                                     ]
    [ -r  | --reverseReads        | pathname for the original reverse fastq reads                                     ]
    [ -o  | --outputDirectory     | directory for output bins                                                         ]
    [ -a  | --assembly            | pathname for the metagenome assembly fasta file to bin                            ]
    [ -b  | --binner              | tool to use for binning; must be one of "metabat", "maxbin", or "concoct"         ]
    [ -l  | --minLength           | minimum contig length to bin; the default is 2000bp                               ]
    [ -e  | --metawrapEnvironment | metaWRAP environment; default "/data/biocore/programs/mamba-envs/metawrap-env"    ]
    [ -c  | --checkm2Environment  | checkM2 environment; default "/data/biocore/programs/mamba-envs/checkm2-anaconda" ]
    [ -h  | --help                | prints an informational message and exits script                                  ]
EOF
  exit;
fi

#### run metaWrap binning tool with specified binner #####

source activate "$mwEnv"

if [[ "$binner" == "maxbin" ]]; then
  metawrap binning \
    -a "$assembly" \
    -o "$outDir" \
    -t 32 -m 256 -l $minLength \
    --maxbin2 --universal \
    "$forwards" "$reverses"
  binDir="$outDir"/maxbin2_bins
  checkDir="$outDir"/maxbin2_checkM2
elif [[ "$binner" == "metabat" ]]; then
  metawrap binning \
    -a "$assembly" \
    -o "$outDir" \
    -t 32 -m 256 -l $minLength \
    --metabat2 \
    "$forwards" "$reverses"
  binDir="$outDir"/metabat2_bins
  checkDir="$outDir"/metabat2_checkM2
elif [[ "$binner" == "concoct" ]]; then
  metawrap binning \
    -a "$assembly" \
    -o "$outDir" \
    -l $minLength \
    -t 32 -m 256 \
    --concoct \
    "$forwards" "$reverses"
  binDir="$outDir"/concoct_bins
  checkDir="$outDir"/concoct_checkM2
else
  echo "no valid binning tool specified"
fi

source deactivate

source activate "$chEnv"

checkm2 predict \
  --threads 32 -x "fa" \
  --input "$binDir" \
  --output-directory "$checkDir"
