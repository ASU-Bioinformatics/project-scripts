#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 0-8:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=128G

module purge
module load mamba/latest

##### Define Variables #####
umask 0007

qiimeDir="$pwd"/qiime2
metadata="$pwd"/metadata.txt
environment="/data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7/"
classifier="/data/biocore/qiime2_classifiers/qiime-2025.7/2024.09.greengenes.backbone.full-length.nb.sklearn-1.4.2.qza"
samplingDepth=10000
minDepth=100
phylogeny="TRUE"
inputStrand="paired"
dada2="single"
help="FALSE"

VALID_ARGS=$(getopt -o q:m:c:n:x:s:d:r:pe:h \
                    --long qiimeDir:,metadata:,categoricals:,numericals:,samplingDepth:,strandedness:,minDepth:,classifier:,phylogeny,environment:,help \
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
    -c | --categoricals)
        echo "The metadata columns to compare are '$2'"
        categoricals=$2
        shift 2
        ;;
    -n | --numericals)
        numericals=$2
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
    -x | --samplingDepth)
        echo "Rarefaction curve sampling will max at '$2'"
        shift 2
        ;;
    -d | --minDepth)
        echo "Rarefaction curve sampling will start at '$2'"
        shift 2
        ;;
    -r | --classifier)
        echo "The classifier to use for taxonomic analysis is '$2'"
        classifier="$2"
        shift 2
        ;;
    -p | --phylogeny)
        echo "Phylogenetic analysis will be excluded '$2'"
        phylogeny="FALSE"
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
    [ -c  |   --categoricals  |   list of column names in the metadata that should be analyzed categorically                                                            ]
    [ -n  |   --numericals    |   list of column names in the metadata that should be analyzed numerically                                                              ]
    [ -s  |   --pairing       |   the pairing strategy of the sequencing. Allowable options are p, s, and ps (paired, single, and paired load with single DADA2)        ]
    [ -x  |   --samplingDepth |   the lowest feature count from the samples table (or from the included samples in that table)                                          ]
    [ -d  |   --minDepth      |   the lowest number of reads to sample for feature analysis curve towards sampling depth; default is 100 reads                          ]
    [ -r  |   --classifier    |   the set of reference sequences used to taxonomically classify the reads; default is the newest version of Greengenes full length      ]
    [ -p  |   --noPhylogeny   |   if selected, this turns off the phylogeny analysis to allow taxonomic analysis to be rerun independently with an alternate classifer  ]
    [ -e  |   --environment   |   location for the Qiime2 environment to activate                                                                                       ]
    [ -h  |   --help          |   prints an informational message and exits script                                                                                      ]
EOF
  exit;
fi

##### Run Analysis #####

source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7/
cd "$qiimeDir"

# Taxonomic analysis with selected classifier
echo $classifier
qiime feature-classifier classify-sklearn \
  --i-classifier $classifier \
  --i-reads rep-seqs-"$inputStrand"-"$dada2".qza \
  --o-classification taxonomy-"$inputStrand"-"$dada2".qza

qiime metadata tabulate \
  --m-input-file taxonomy-"$inputStrand"-"$dada2".qza \
  --o-visualization taxonomy-"$inputStrand"-"$dada2".qzv

qiime taxa barplot \
  --i-table table-"$inputStrand"-"$dada2".qza \
  --i-taxonomy taxonomy-"$inputStrand"-"$dada2".qza \
  --m-metadata-file "$metadata" \
  --o-visualization taxa-bar-plots-"$inputStrand"-"$dada2".qzv

for i in 1 2 3 4 5 6 7
do
  qiime taxa collapse \
    --i-table table-"$inputStrand"-"$dada2".qza \
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

#phylogeny analysis (classifier neutral)
if [ $phylogeny == "TRUE" ];
then
  qiime phylogeny align-to-tree-mafft-fasttree \
    --i-sequences rep-seqs-"$inputStrand"-"$dada2".qza \
    --o-alignment aligned-rep-seqs-"$inputStrand"-"$dada2".qza \
    --o-masked-alignment masked-aligned-rep-seqs-"$inputStrand"-"$dada2".qza \
    --o-tree unrooted-tree-"$inputStrand"-"$dada2".qza \
    --o-rooted-tree rooted-tree-"$inputStrand"-"$dada2".qza

  qiime diversity core-metrics-phylogenetic \
    --i-phylogeny rooted-tree-"$inputStrand"-"$dada2".qza \
    --i-table table-"$inputStrand"-"$dada2".qza \
    --p-sampling-depth $samplingDepth \
    --m-metadata-file "$metadata" \
    --output-dir core-metrics-results-"$inputStrand"-"$dada2"

  qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results-"$inputStrand"-"$dada2"/faith_pd_vector.qza \
    --m-metadata-file "$metadata" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/faith-pd-group-significance.qzv

  qiime diversity alpha-group-significance \
    --i-alpha-diversity core-metrics-results-"$inputStrand"-"$dada2"/evenness_vector.qza \
    --m-metadata-file "$metadata" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/evenness-group-significance.qzv

#beta group significance, for each categorical column
for j in $categoricals;
do

  (echo "$j"

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/unweighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/weighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"/bray_curtis_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/bray_curtis_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"/jaccard_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/jaccard_"$j"_significance.qzv \
    --p-pairwise

  ) &

done;
wait

for k in $numericals;
do

  (echo "$k"

  qiime emperor plot \
    --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"/unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file "$metadata" \
    --p-custom-axes "$k" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/unweighted_unifrac_emperor-"$k".qzv \

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"/bray_curtis_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/bray-curtis-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"/jaccard_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/jaccard-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"/weighted_unifrac-"$k".qzv

  ) &

done;
wait


# Alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table table-"$inputStrand"-"$dada2".qza \
  --i-phylogeny rooted-tree-"$inputStrand"-"$dada2".qza \
  --p-min-depth $minDepth \
  --p-max-depth $samplingDepth \
  --p-steps 200 \
  --p-iterations 10 \
  --m-metadata-file "$metadata" \
  --o-visualization alpha-rarefaction-"$inputStrand"-"$dada2".qzv
fi
