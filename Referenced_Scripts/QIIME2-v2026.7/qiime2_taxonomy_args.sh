#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.q2tax.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.q2tax.err                   # STDERR (%j = JobId)
#SBATCH -t 0-4:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=96G

module purge
module load mamba/latest

##### Define Variables #####
umask 0007

qiimeDir="$pwd"/qiime2
metadata="$pwd"/metadata.txt
environment="/data/biocore/programs/mamba-envs/qiime2-2026.7/"
classifier="/data/biocore/qiime2_classifiers/qiime-2026.7/SILVA_144_SSURef_NR99_uniform_classifier_full-length.qza"
inputStrand="paired"
dada2="single"
help="FALSE"

VALID_ARGS=$(getopt -o q:m:s:r:e:h \
                    --long qiimeDir:,metadata:,strandedness:,classifier:,environment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -q | --qiimeDir)
        echo "Output Qiime2 files will be written to '$2'"
        qiimeDir="$2"
        shift 2
        ;;
    -m | --metadata)
        echo "The sample metadata file to use is '$2'"
        metadata="$2"
        shift 2
        ;;
    -s | --strandedness)
        if [ "$2" == "p" ];
          then
            echo "The read type is paired-end"
            inputStrand="paired"
            dada2="paired"
          elif [ "$2" == "s" ];
          then
            echo "The read type is single-end"
            inputStrand="single"
            dada2="single"
          elif [ "$2" == "ps" ];
          then
            echo "The read type is paired-end but DADA2 will be run with single-end parameters"
            inputStrand="paired"
            dada2="single"
        fi
        shift 2
        ;;
    -r | --classifier)
        echo "The classifier to use for taxonomic analysis is '$2'"
        classifier="$2"
        class=$(basename "$classifier")
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
  This script runs the first part of our standard Qiime2 analysis (loading the data, denoising,
  and summarizing statistics by metadata information).

  usage: sbatch qiime2_pt2_args.sh
            -q /path/to/qiime-output -m /path/to/metadata.txt
            -c "list of categorical columns" -n "list of numerical columns"
            -s ps -x 10000 -d 100
            -r 2024.09.greengenes.backbone.full-length.nb.sklearn-1.4.2.qza
            -e /path/to/conda/environment (-p) (-h)

  options:
    [ -q  |   --qiimeDir      |   directory for Qiime2 output files (will be created if it doesn't already exist; previous files will be overwritten)                  ]
    [ -m  |   --metadata      |   text file containing metadata information in Qiime2-compatible format                                                                ]
    [ -s  |   --pairing       |   the pairing strategy of the sequencing. Allowable options are p, s, and ps (paired, single, and paired load with single DADA2)        ]
    [ -r  |   --classifier    |   the set of reference sequences used to taxonomically classify the reads; default is the newest version of Greengenes full length      ]
    [ -e  |   --environment   |   location for the Qiime2 environment to activate                                                                                       ]
    [ -h  |   --help          |   prints an informational message and exits script                                                                                      ]
EOF
  exit;
fi

##### Run Analysis #####

source activate "$environment"

mkdir -p "$qiimeDir"/"$class"-taxonomy
cd "$qiimeDir"/"$class"-taxonomy

# Taxonomic analysis with selected classifier
echo $classifier
qiime feature-classifier classify-sklearn \
  --i-classifier $classifier \
  --i-reads "$qiimeDir"/rep-seqs-"$inputStrand"-"$dada2".qza \
  --o-classification taxonomy-"$inputStrand"-"$dada2".qza

qiime metadata tabulate \
  --m-input-file taxonomy-"$inputStrand"-"$dada2".qza \
  --o-visualization taxonomy-"$inputStrand"-"$dada2".qzv

qiime taxa barplot \
  --i-table "$qiimeDir"/table-"$inputStrand"-"$dada2".qza \
  --i-taxonomy taxonomy-"$inputStrand"-"$dada2".qza \
  --m-metadata-file "$metadata" \
  --o-visualization taxa-bar-plots-"$inputStrand"-"$dada2".qzv

for i in 1 2 3 4 5 6 7
do
  qiime taxa collapse \
    --i-table "$qiimeDir"/table-"$inputStrand"-"$dada2".qza \
    --i-taxonomy taxonomy-"$inputStrand"-"$dada2".qza \
    --o-collapsed-table level"$i"_table-"$inputStrand"-"$dada2".qza \
    --p-level $i

  qiime feature-table relative-frequency \
    --i-table level"$i"_table-"$inputStrand"-"$dada2".qza \
    --o-relative-frequency-table rel-level"$i"_table-"$inputStrand"-"$dada2".qza

  mkdir rel-table"$i"-"$inputStrand"-"$dada2"
  qiime tools export \
    --input-path rel-level"$i"_table-"$inputStrand"-"$dada2".qza \
    --output-path rel-table"$i"-"$inputStrand"-"$dada2"
  #this step will probably produce a python error but it succeeds
  cd rel-table"$i"-"$inputStrand"-"$dada2"
  biom convert -i feature-table.biom -o rel-level"$i"-table-"$inputStrand"-"$dada2".tsv --to-tsv
  cd ../
done
