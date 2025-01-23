#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.%a.align.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.%a.align.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 1-0:00
#SBATCH -c 6
#SBATCH --mem=32G

module purge

paired="TRUE"
help="FALSE"
itmdir=$(pwd)

VALID_ARGS=$(getopt -o s:f:i:p:r:o:ta:h \
                    --long sampleID:,fastqDir:,intermediateDir:,starParams:,refDir:,outputDir:,readTypeSingle,array:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -s | --sampleID)
        echo "Aligning sample '$2'"
        sid="$2"
        shift 2
        ;;
    -f | --fastqDir)
        echo "The sample can be found in the directory '$2'"
        fastqDir="$2"
        shift 2
        ;;
    -i | --intermediateDir)
        echo "FIFO server drive for intermediate alignment is '$2'"
        itmdir="$2"
        shift 2
        ;;
    -p | --starParams)
        echo "STAR alignment parameter text file is '$2'"
        params="$2"
        shift 2
        ;;
    -r | --refDir)
        echo "Reference files and STAR indexes are in '$2'"
        rdir="$2"
        shift 2
        ;;
    -o | --outputDir)
        echo "Results will be output to the directory '$2'"
        outputDir="$2"
        shift 2
        ;;
    -t | --readTypeSingle)
        echo "Samples are from single-end read sequencing"
        paired="FALSE"
        shift
        ;;
    -a | --array)
        array=$2
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
  This script is used to quantify the abundance of each transcript
  from the reference GTF in a single STAR alignment file.

  usage: bash A1.cpu-align.sh -s id -f fdir -i idir -p txt -r rdir -o odir (-t) (-d)

  options:
    [ -s  |   --sampleID        |   string before the first underscore in the fastq.gz file name          ]
    [ -f  |   --fastqDir        |   directory in which samples can be found                               ]
    [ -i  |   --intermediateDir |   FIFO-enabled server for STAR alignment                                ]
    [ -p  |   --starParams      |   text file containing STAR parameters (see examples in scripts folder) ]
    [ -r  |   --refDir          |   directory containing the reference genome, GTF file, and STAR indexes ]
    [ -o  |   --outputDir       |   location for output alignment files (eg, sorted bam files)            ]
    [ -t  |   --readTypeSingle  |   optional; specifies that input reads are single-end (unidirectional)  ]
    [ -a  |   --array           |   DO NOT USE; this enables connection to workflow wrapper scripts       ]
    [ -h  |   --help            |   prints an informational message and exists                            ]

EOF
  exit;
fi

if [[ -v array ]]
then
  samples=($array)
  sid=${samples[$SLURM_ARRAY_TASK_ID]}
fi

mkdir -p "$itmdir"
mkdir -p "$outputDir"

echo $sid

cp "$fastqDir"/"$sid"_S*.fastq.gz "$itmdir"

cd "$itmdir"

ulimit -v 30000000000

if [ "$paired" == "TRUE" ]
then
  ls -lh "$sid"_S*_R1_*.fastq.gz;
  ls -lh "$sid"_S*_R2_*.fastq.gz;
  STAR \
    --genomeDir "$rdir" \
    --readFilesCommand gunzip -c \
    --readFilesIn "$sid"_S*_R1_*.fastq.gz "$sid"_S*_R2_*.fastq.gz \
    --outFileNamePrefix "$sid"_STAR \
    --limitBAMsortRAM 26021018803 \
    --parametersFiles "$params"
else
  ls -lh $i;
  STAR \
    --genomeDir "$rdir" \
    --readFilesCommand gunzip -c \
    --readFilesIn "$sid"_S*_R1_*.fastq.gz \
    --outFileNamePrefix "$sid"_STAR \
    --limitBAMsortRAM 26021018803 \
    --parametersFiles "$params"
fi

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

rm "$itmdir"/"$sid"_S*.fastq.gz

mv "$sid"_S* "$outputDir"
