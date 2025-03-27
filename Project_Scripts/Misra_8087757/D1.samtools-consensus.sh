#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.consensus.out
#SBATCH -e slurm.%A.consensus.err
#SBATCH -t 0-8:00
#SBATCH -c 1
#SBATCH --mem=64G

VALID_ARGS=$(getopt -o i:o:h \
                    --long inBam:,outFasta:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inBam)
        inBam="$2"
        shift 2
        ;;
    -o | --outFasta)
        outFasta="$2"
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
  This script builds a reference-type genome sequence based off of an aligned bam file.
  It includes both insertions and deletions, marking insertions so that the coordinates
  still match the coordinates in the reference genome.

  usage: bash D1.samtools-consensus.sh \
                -i /path/to/inBam \
                -o /path/to/outFasta (-h)

  options:
    [ -i  |   --inBam      |   bam file of reads aligned to a reference genome  ]
    [ -o  |   --outFasta   |   name of the output consensus fasta sequence      ]
    [ -h  |   --help       |   prints an informational message and exits script ]
EOF
  exit;
fi

module load samtools-1.21-gcc-12.1.0

samtools consensus \
  -a --show-del --show-ins --mark-ins \
  -o "$outFasta" "$inBam"
