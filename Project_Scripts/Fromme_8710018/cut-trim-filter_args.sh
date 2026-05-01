#!/bin/bash

#SBATCH -p general
#SBATCH -q public                     # sol QOS
###SBATCH -q grp_kawoodbu           # phx QOS
#SBATCH -o slurm.%A.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-12:00
#SBATCH --mem=32G

# this script is used for removing adapters; no quality control is performed
# the sbatch parameters are set for phx by default; if running on Sol change sbatch QOS

#the script is designed to take any input fastq.gz files where the sample ID does not contain and underscore and is separated from the remainder of
# the file name by an underscore, as long as a R1 or R2 can be found in between the sample ID and extension.
# the output files, however, will be formatted with the traditional Illumina/Casava naming system.
umask 0007

inputDir="$pwd"/fastq
cutadaptDir="$pwd"/fastq-cutadapt
adapters="/data/gencore/databases/trimmomatic/all.fa"
environment="/data/biocore/programs/mamba-envs/cutadapt/"

VALID_ARGS=$(getopt -o i:c:a:m:l:e:h \
                    --long inputDir:,cutadaptDir:,adapters:,cutParams:,cutLen:,environment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Input fastq files will be read from '$2'"
        inputDir="$2"
        shift 2
        ;;
    -c | --cutadaptDir)
        echo "Intermediate fastq files output by cutadapt only will be written to '$2'"
        cutadaptDir="$2"
        shift 2
        ;;
    -a | --adapters)
        echo "The fastq file of adapters to trim is '$2'"
        shift 2
        ;;
    -e | --environment)
        echo "The conda environment to use is '$2'"
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
  This script removes adapter sequences and filters for quality in chloroplast reads. Honestly, it could be used
  for a variety of read types beyond that. The cutadapt parameters can't be customized at this time, and are -m 30 -q 20

  usage: sbatch cut-trim-filter_args.sh
            -i /path/to/fastq-input
            -c /path/to/cutadapt-output
            -a /path/to/adapters.fa
            -e /path/to/conda/environment (-h)

  options:
    [ -i  |   --inputDir     |   directory containing fastq.gz files, where the sample ID is the first field before and underscore  ]
                                   and the read designation (R1 or R2) is in between the sample ID and extension                    ]
    [ -c  |   --cutadaptDir  |   directory for files output by cutadapt, before trimmomatic quality filtering                       ]
    [ -a  |   --adapters     |   the file containing fastq sequences for adapters to search for and trim                            ]
    [ -e  |   --environment  |   location for the cutadapt environment; default is "/data/biocore/programs/mamba-envs/cutadapt/"    ]
    [ -h  |   --help         |   prints an informational message and exits script                                                   ]
EOF
  exit;
fi

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/cutadapt/

mkdir -p $cutadaptDir
cd $inputDir

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  cutadapt -a file:"$adapters" -g file:"$adapters" \
         -A file:"$adapters" -G file:"$adapters" \
         -m 30 -q 20 \
         -o "$i"_SCT_L001_R1_001.fastq.gz \
         -p "$i"_SCT_L001_R2_001.fastq.gz \
         "$i"_S*R1*.fastq.gz "$i"_S*R2*.fastq.gz
done

mv *SCT* $cutadaptDir
