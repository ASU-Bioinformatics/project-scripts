#!/bin/bash

#SBATCH -p public                    # sol partition
####SBATCH -p general                # phx partition
#SBATCH -q public                    # sol QOS
####SBATCH -q grp_kawoodbu           # phx QOS
#SBATCH -o slurm.%.cut.out
#SBATCH -e slurm.%.cut.err
#SBATCH -t 2-0:00
#SBATCH --mem=64G
#SBATCH -c 4

# this is the script is used for removing adapters and trimming for quality
# the sbatch parameters are set for Sol by default; if running on phx change sbatch QOS and trimmomatic module

#the script is designed to take any input fastq.gz files where the sample ID does not contain an underscore and is separated from the remainder of
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
cut=TRUE
trim=TRUE

VALID_ARGS=$(getopt -o i:c:a:e:h \
                    --long inputDir:,cutadaptDir:,adapters:,environment:,help \
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
        echo "Fastq files output by cutadapt will be written to '$2'"
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
  This script removes adapter sequences and filters for quality in 18S, 16S, and ITS reads. Honestly, it could be used
  for a variety of read types beyond that. The cutadapt parameters can't be customized at this time, and are -m 30 -q 10

  usage: sbatch cut-trim-filter_args.sh
            -i /path/to/fastq-input
            -c /path/to/cutadapt-output
            -a /path/to/adapters.fa
            -e /path/to/conda/environment (-h)

  options:
    [ -i  |   --inputDir      |   directory containing fastq.gz files, where the sample ID is the first field before an underscore   ]
                                   and the read designation (R1 or R2) is in between the sample ID and extension                     ]
    [ -c  |   --cutadaptDir   |   directory for files output by cutadapt, before trimmomatic quality filtering                       ]
    [ -a  |   --adapters      |   the file containing fastq sequences for adapters to search for and trim                            ]
    [ -e  |   --environment   |   location for the cutadapt environment; default is "/data/biocore/programs/mamba-envs/cutadapt/"    ]
    [ -h  |   --help          |   prints an informational message and exits script                                                   ]
EOF
  exit;
fi

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/cutadapt/

mkdir -p $cutadaptDir
cd $inputDir

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq )
do
  ( echo "$i"
  cutadapt -a file:"$adapters" -A file:"$adapters" \
         -m 30 -q 15 --cores=4 \
         -o "$i"_SCT_L001_R1_001.fastq.gz \
         -p "$i"_SCT_L001_R2_001.fastq.gz \
         "$i"_S*R1*.fastq.gz "$i"_S*R2*.fastq.gz
       ) &
done;
wait

mv *SCT_L001_R1_001.fastq.gz $cutadaptDir
mv *SCT_L001_R2_001.fastq.gz $cutadaptDir

cd $cutadaptDir
