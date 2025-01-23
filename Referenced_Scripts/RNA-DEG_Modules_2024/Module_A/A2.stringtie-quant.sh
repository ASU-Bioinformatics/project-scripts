#!/bin/bash

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%A.%a.out                # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.%a.err                # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-1:00                         # estimated time needed
#SBATCH --mem=8G

module purge
module load mamba/latest

source activate /data/biocore/programs/mamba-envs/biocore-rna

PATH=/data/biocore/programs/mamba-envs/biocore-rna/bin/:$PATH

help="FALSE"

VALID_ARGS=$(getopt -o s:q:g:a:l:h \
                    --long sampleID:,quantDir:,refGTF:,alignmentDir:,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -s | --sampleID)
        echo "Quantifying alignents from sample '$2'"
        sid="$2"
        shift 2
        ;;
    -q | --quantDir)
        echo "Output stringtie quant files to '$2'"
        qdir="$2"
        shift 2
        ;;
    -g | --refGTF)
        echo "The reference GTF file is '$2'"
        gtf="$2"
        shift 2
        ;;
    -a | --alignmentDir)
        echo "STAR alignment bam files are found in '$2'"
        adir="$2"
        shift 2
        ;;
    -l | --list)
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

  usage: bash A2.stringtie-quant.sh -s id -a adir -q qdir -g gtf

  options:
    [ -s  |   --sampleID      |   string before the first underscore in the bam file name         ]
    [ -a  |   --alignmentDir  |   directory containing STAR alignment bam files                   ]
    [ -q  |   --quantDir      |   directory for stringtie output                                  ]
    [ -g  |   --refGTF        |   full pathname to the GTF file used for alignment                ]
    [ -l  |   --list          |   DO NOT USE; this enables connection to workflow wrapper scripts ]

EOF
exit;
fi

mkdir -p "$qdir"

if [[ -v array ]]
then
  samples=($array)
  sid=${samples[$SLURM_ARRAY_TASK_ID]}
fi

#stringtie quantification
if [ -d "$qdir"/"$sid"/tmp_* ]; then
  echo "stringtie temp directory detected"
  rm -r "$qdir"/"$sid"
fi

stringtie \
  -e -B -p 8 -G "$gtf" -v \
  -A "$qdir"/"$sid"/"$sid".stringtie.abund.txt \
  -o "$qdir"/"$sid"/"$sid".stringtie.gtf \
  "$adir"/"$sid"_STARAligned.sortedByCoord.out.bam
