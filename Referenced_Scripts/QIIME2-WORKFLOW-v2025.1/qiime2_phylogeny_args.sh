#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.q2phylo.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.q2phylo.err                   # STDERR (%j = JobId)
#SBATCH -t 0-8:00                         
#SBATCH --mem=32G

module purge
module load mamba/latest

##### Define Variables #####
umask 0007

qiimeDir="$pwd"/qiime2
metadata="$pwd"/metadata.txt
environment="/data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7/"
samplingDepth=10000
minDepth=100
inputStrand="paired"
dada2="paired"
help="FALSE"

VALID_ARGS=$(getopt -o q:m:c:n:x:s:d:e:h \
                    --long qiimeDir:,metadata:,categoricals:,numericals:,samplingDepth:,strandedness:,minDepth:,environment:,help \
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
        samplingDepth="$2"
        shift 2
        ;;
    -d | --minDepth)
        echo "Rarefaction curve sampling will start at '$2'"
        minDepth="$2"
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
  This script runs the metadata-and-phylogeny-based analysis of microbial samples,
  following dada2 denoising.

  usage: sbatch qiime2_pt2_phylogeny_args.sh
            -q /path/to/qiime-output -m /path/to/metadata.txt
            -c "list of categorical columns" -n "list of numerical columns"
            -s ps -x 10000 -d 100
            -e /path/to/conda/environment (-p) (-h)

  options:
    [ -q  |   --qiimeDir      |   directory for Qiime2 output files (will be created if it doesn't already exist; previous files will be overwritten)                  ]
    [ -m  |   --metadata      |   text file containing metadata information in Qiime2-compatible format                                                                ]
    [ -c  |   --categoricals  |   list of column names in the metadata that should be analyzed categorically                                                            ]
    [ -n  |   --numericals    |   list of column names in the metadata that should be analyzed numerically                                                              ]
    [ -s  |   --pairing       |   the pairing strategy of the sequencing. Allowable options are p, s, and ps (paired, single, and paired load with single DADA2)        ]
    [ -x  |   --samplingDepth |   the lowest feature count from the samples table (or from the included samples in that table)                                          ]
    [ -d  |   --minDepth      |   the lowest number of reads to sample for feature analysis curve towards sampling depth; default is 100 reads                          ]
    [ -e  |   --environment   |   location for the Qiime2 environment to activate                                                                                       ]
    [ -h  |   --help          |   prints an informational message and exits script                                                                                      ]
EOF
  exit;
fi

##### Run Analysis #####

source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7/
cd "$qiimeDir"

#phylogeny analysis (classifier neutral)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs-"$inputStrand"-"$dada2".qza \
  --o-alignment aligned-rep-seqs-"$inputStrand"-"$dada2"-"$samplingDepth".qza \
  --o-masked-alignment masked-aligned-rep-seqs-"$inputStrand"-"$dada2"-"$samplingDepth".qza \
  --o-tree unrooted-tree-"$inputStrand"-"$dada2"-"$samplingDepth".qza \
  --o-rooted-tree rooted-tree-"$inputStrand"-"$dada2"-"$samplingDepth".qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree-"$inputStrand"-"$dada2"-"$samplingDepth".qza \
  --i-table table-"$inputStrand"-"$dada2".qza \
  --p-sampling-depth $samplingDepth \
  --m-metadata-file "$metadata" \
  --output-dir core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/faith_pd_vector.qza \
  --m-metadata-file "$metadata" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/evenness_vector.qza \
  --m-metadata-file "$metadata" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/evenness-group-significance.qzv

#beta group significance, for each categorical column
for j in $categoricals;
do

  (echo "$j"

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/unweighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/weighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/bray_curtis_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/bray_curtis_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/jaccard_distance_matrix.qza \
    --m-metadata-file "$metadata" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/jaccard_"$j"_significance.qzv \
    --p-pairwise

  ) &

done;
wait

for k in $numericals;
do

  (echo "$k"

  qiime emperor plot \
    --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file "$metadata" \
    --p-custom-axes "$k" \
    --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/unweighted_unifrac_emperor-"$k".qzv \

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/bray_curtis_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/bray-curtis-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/jaccard_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/jaccard-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file "$metadata" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results-"$inputStrand"-"$dada2"-"$samplingDepth"/weighted_unifrac-"$k".qzv

  ) &

done;
wait


# Alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table table-"$inputStrand"-"$dada2".qza \
  --i-phylogeny rooted-tree-"$inputStrand"-"$dada2"-"$samplingDepth".qza \
  --p-min-depth $minDepth \
  --p-max-depth $samplingDepth \
  --p-steps 200 \
  --p-iterations 10 \
  --m-metadata-file "$metadata" \
  --o-visualization alpha-rarefaction-"$inputStrand"-"$dada2"-"$samplingDepth".qzv
