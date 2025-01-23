#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.merge.out                # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.merge.err                # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-4:00                         # estimated time needed

module purge
module load mamba/latest

source activate /data/biocore/programs/mamba-envs/biocore-rna
PATH=/data/biocore/programs/mamba-envs/biocore-rna/bin/:$PATH

help="FALSE"

VALID_ARGS=$(getopt -o s:q:a:h \
                    --long scriptDir:,quantDir:,alignmentDir:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -s | --scriptDir)
        echo "prepDE.py and stringtie_expression_matrix.pl scripts are located in '$2'"
        sdir="$2"
        shift 2
        ;;
    -q | --quantDir)
        echo "Stringtie input folders are located in '$2'"
        qdir="$2"
        shift 2
        ;;
    -a | --alignmentDir)
        echo "STAR alignments to summarize are located in '$2'"
        adir="$2"
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
  This script is used to merge stringtie quantifications from a sample set,
  creating (among other things) a single gene count matrix.

  The source directory should contain prepDE.py, stringtie_expression_matrix.pl,
  clean-counts.py, and star-summary.py.

  usage: bash merge-quant.sh -s scriptDir -q quantDir

  options:
    [ -s  |   --scriptDir     |   directory containing the auxiliary perl and python scripts  ]
    [ -q  |   --quantDir      |   directory containing the sample folders from stringtie      ]
    [ -a  |   --alignmentDir  |   directory containing the STAR alignment files to summarize  ]

EOF
exit;
fi

cd "$qdir"

python2 "$sdir"/A4.prepDE.py -i . #if python3 is saved as python in your path, specify python2 here

for i in $(find "$qdir" -mindepth 1 -type d | while read F; do basename $F; done)
do
  cd "$qdir"/"$i"
  cp "$i".stringtie.abund.txt gene_abundances.tsv
  cp "$i".stringtie.gtf transcripts.gtf
done;

cd "$qdir"

#store list of directories in a variable for passing to stringtie_expression_matrix.pl
dirs=$(find ./ -type d | paste -d, -s | sed 's/\.//g' | sed 's/\///g' | sed 's/^.//')
echo 'Building an expression matrix for the following directories:'
echo $dirs

export PERL5LIB=/data/biocore/programs/miniconda3/bin/perl

perl "$sdir"/A5.stringtie_expression_matrix.pl \
  --expression_metric=TPM \
  --result_dirs=$dirs \
  --transcript_matrix_file=transcript_tpm_all_samples.tsv \
  --gene_matrix_file=gene_tpm_all_samples.tsv

perl "$sdir"/A5.stringtie_expression_matrix.pl \
  --expression_metric=FPKM \
  --result_dirs=$dirs \
  --transcript_matrix_file=transcript_fpkm_all_samples.tsv \
  --gene_matrix_file=gene_fpkm_all_samples.tsv

python "$sdir"/A6.clean-counts.py gene_count_matrix.csv

cd "$adir"

python "$sdir"/A7.star-summary.py "$adir"
