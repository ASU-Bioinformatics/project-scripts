#!/bin/bash

#### fastqc file generation #####

#SBATCH -p general
#SBATCH -q grp_kawoodbu #phx
#SBATCH -t 1-0:00               # estimated time needed
#SBATCH --mem=256G
#SBATCH -c 2

umask 0007

dir="TRUE"
environment="/data/biocore/programs/mamba-envs/multiqc.v1.20/"

VALID_ARGS=$(getopt -o f:i:e:h \
                    --long fastqDir:,inputFiles:,environment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --fastqDir)
        echo "fastqc will be run on all fastq.gz files in the directory '$2'"
        fastqDir="$2"
        shift 2
        ;;
    -i | --inputFiles)
        echo "fastqc will be run on all fastq.gz files listed in the file '$2'"
        inputFiles="$2"
        dir="FALSE"
        shift 2
        ;;
    -e | --environment)
        echo "The conda environment to use for multiqc is '$2'"
        environment="$2"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
        break
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
  This script generates fastqc files for all files in a directory and summarizes them with multiqc.
  Alternately, fastqc and multiqc can be run only for files specified in a manifest file. In this case,
  the fastqDir parameter should still be included and will serve as the output directory.
  The conda environment can be specified as well as the fastq directory.

  usage: sbatch fastqc_generation_args_phx.sh
            -f /path/to/fastq-input
            -i /path/to/input-files.txt
            -e /path/to/conda/environment (-h)

  options:
    [ -f  |   --fastqDir     |   directory containing fastq.gz files, where the sample ID is the first field before and underscore    ]
    [ -i  |   --inputFiles   |   a file containing a list of all fastq.gz files to qc, one file name per line                         ]
    [ -e  |   --environment  |   location for the multiqc environment; default is "/data/biocore/programs/mamba-envs/multiqc.v1.20/"  ]
    [ -h  |   --help         |   prints an informational message and exits script                                                     ]
EOF
  exit;
fi

module purge
module load fastqc-0.12.1-v5 #phx

##### go to directory where fastq files are located, which is also the output folder #####

cd "$fastqDir"
mkdir -p fastq
mkdir -p qc

if [ "$dir" == "FALSE" ]; then
  fastqc -t 256 $(cat "$inputFiles")
  mv $(cat "$inputFiles") fastq
else
  fastqc -t 256 *.fastq.gz
  mv *fastq.gz fastq
fi

mv *fastqc* qc

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/multiqc.v1.20/

multiqc "$fastqDir/qc"

chmod -R g+w *
