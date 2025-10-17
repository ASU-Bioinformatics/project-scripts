#!/bin/bash

#SBATCH -p general
####SBATCH -q public              # sol QOS
#SBATCH -q grp_kawoodbu           # phx QOS
#SBATCH -o slurm.%A.cut-trim.out
#SBATCH -e slurm.%A.cut-trim.err
#SBATCH -t 2-0:00
#SBATCH --mem=64G

# this is the script is used for removing adapters and trimming for quality
# the sbatch parameters are set for phx by default; if running on Sol change sbatch QOS and trimmomatic module

#the script is designed to take any input fastq.gz files where the sample ID does not contain and underscore and is separated from the remainder of
# the file name by an underscore, as long as a R1 or R2 can be found in between the sample ID and extension.
# the output files, however, will be formatted with the traditional Illumina/Casava naming system.
umask 0007

inputDir="$pwd"/fastq
cutadaptDir="$pwd"/fastq-cutadapt
pairedDir="$pwd"/fastq-cut-paired
unpairedDir="$pwd"/fastq-cut-unpaired
adapters="/data/gencore/databases/trimmomatic/all.fa"
cropLength=300
minLength=50
environment="/data/biocore/programs/mamba-envs/cutadapt/"

VALID_ARGS=$(getopt -o i:c:p:u:a:t:m:l:e:h \
                    --long inputDir:,cutadaptDir:,pairedDir:,unpairedDir:,adapters:,trimParams:,cutParams:,cutLen:,environment:,help \
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
    -p | --pairedDir)
        echo "The paired fastq files output from trimmomatic following cutadapt will be written to '$2'"
        pairedDir="$2"
        shift 2
        ;;
    -u | --unpairedDir)
        echo "The unpaired fastq files output from trimmomatic following cutadapt will be written to '$2'"
        unpairedDir="$2"
        shift 2
        ;;
    -a | --adapters)
        echo "The fastq file of adapters to trim is '$2'"
        shift 2
        ;;
    -t | --trimParams)
        trimParams="$2"
        echo "The trimmomatic parameters are '$2'"
        shift 2
        ;;
    -m | --cutParams)
        cutParams="$2"
        echo "The cutadapt parameters are '$2'"
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
  This script removes adapter sequences and filters for quality in 18S, 16S, and ITS reads. Honestly, it could be used
  for a variety of read types beyond that. The cutadapt parameters can't be customized at this time, and are -m 30 -q 20

  usage: sbatch cut-trim-filter_args.sh
            -i /path/to/fastq-input
            -c /path/to/cutadapt-output
            -p /path/to/paired-trimmomatic-output
            -u /path/to/unpaired-trimmomatic-output
            -a /path/to/adapters.fa
            -t "trimmomatic parameters"
            -e /path/to/conda/environment (-h)

  options:
    [ -i  |   --inputDir     |   directory containing fastq.gz files, where the sample ID is the first field before and underscore  ]
                                   and the read designation (R1 or R2) is in between the sample ID and extension                    ]
    [ -c  |   --cutadaptDir  |   directory for files output by cutadapt, before trimmomatic quality filtering                       ]
    [ -p  |   --pairedDir    |   directory for filtered paired fastq files output by trimmomatic                                    ]
    [ -u  |   --unpairedDir  |   directory for filtered unpaired fastq files output by trimmomatic                                  ]
    [ -a  |   --adapters     |   the file containing fastq sequences for adapters to search for and trim                            ]
    [ -l  |   --cropLength   |   the maximum length to crop reads during Trimmomatic filtering                                      ]
    [ -l  |   --cropLength   |   the minimum read length to retain during Trimmomatic filtering                                     ]
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

mv *SCT_L001_R2_001.fastq.gz $cutadaptDir

cd $cutadaptDir

module load trimmomatic-0.39-fs #phx
###module load trimmomatic-0.39-gcc-12.1.0 #sol

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  trimmomatic PE "$i"_SCT_L001_R1_001.fastq.gz "$i"_SCT_L001_R2_001.fastq.gz \
            "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
            "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
            CROP:$cropLength SLIDINGWINDOW:4:15 MINLEN:$minLength
done

mkdir -p $pairedDir
mkdir -p $unpairedDir

mv *SQP* $pairedDir
mv *SUN* $unpairedDir
