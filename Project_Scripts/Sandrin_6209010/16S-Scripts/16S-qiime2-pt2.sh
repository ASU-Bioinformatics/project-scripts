#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -A kawoodbu
#SBATCH -o slurm.%j.out                 # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                 # STDERR (%j = JobId)
#SBATCH --mail-type=ALL                 # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=kawoodbu@asu.edu  # send-to address
#SBATCH -t 0-8:00                       # estimated time needed - silva classifier takes a long time
#SBATCH --mem=32G

module purge

##### Define Variables #####

METADATA="/data/gencore/analysis_projects/completed_projects/5920921_Sandrin_16S/metadata.txt"
DIRECTORY="/data/gencore/analysis_projects/5920921_Sandrin_16S/qiime2"
SAMPLING_DEPTH=10759 #typically, lowest feature count from samples, found in table.qzv
MIN_DEPTH=100 #something lower, to get a good curve from min_depth to sampling_depth

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2022.2/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2022.2/gg-13-8-99-515-806-nb-classifier.qza"
#classifiers[ez]="/data/biocore/qiime2_classifiers/ezbiocloud_2021_03_17_nb_classifier.qza"
#classifiers[18]="/data/biocore/qiime2_classifiers/sv18s-rep-set-99-nb-classifier.qza"
#classifiers[it]="/data/biocore/qiime2_classifiers/unITs-rep-set-99-nb-classifier.qza"

#put all categorical columns here for unifrac plots by category
categories="DiseaseState
            Sex
            Race
            Diagnosis"

numericals="Age"

##### Run Analysis #####

source activate /data/biocore/programs/conda-envs/qiime2-2022.2

cd "$DIRECTORY"

# Taxonomic analysis with Silva, if selected
for i in ${!classifiers[@]};
do

(echo "$i"

qiime feature-classifier classify-sklearn \
  --i-classifier ${classifiers["$i"]} \
  --i-reads rep-seqs.qza \
  --o-classification "$i"-taxonomy.qza

qiime metadata tabulate \
  --m-input-file "$i"-taxonomy.qza \
  --o-visualization "$i"-taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy "$i"-taxonomy.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization "$i"-test-taxa-bar-plots.qzv

qiime taxa collapse \
    --i-table table.qza \
    --i-taxonomy "$i"-taxonomy.qza \
  --o-collapsed-table "$i"-level7_table.qza \
    --p-level 7

qiime feature-table relative-frequency \
--i-table "$i"-level7_table.qza \
--o-relative-frequency-table "$i"-rel-level7_table.qza

mkdir "$i"-rel-table7
qiime tools export \
--input-path "$i"-rel-level7_table.qza \
--output-path "$i"-rel-table7
#this step will probably produce a python error but it succeeds
cd "$i"-rel-table7
biom convert -i feature-table.biom -o "$i"-rel-level7-table.tsv --to-tsv
cd ..) &

done;
wait

#phylogeny analysis (classifier neutral)
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences rep-seqs.qza \
  --o-alignment aligned-rep-seqs.qza \
  --o-masked-alignment masked-aligned-rep-seqs.qza \
  --o-tree unrooted-tree.qza \
  --o-rooted-tree rooted-tree.qza

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file "$METADATA" \
  --output-dir core-metrics-results

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/faith_pd_vector.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization core-metrics-results/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity core-metrics-results/evenness_vector.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization core-metrics-results/evenness-group-significance.qzv

#beta group significance, for each categorical column
for j in $categories;
do

  (echo "$j"

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results/unweighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results/weighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/bray_curtis_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results/bray_curtis_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix core-metrics-results/jaccard_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization core-metrics-results/jaccard_"$j"_significance.qzv \
    --p-pairwise

  ) &

done;
wait

for k in $numericals;
do

  (echo "$k"

  qiime emperor plot \
    --i-pcoa core-metrics-results/unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file "$METADATA" \
    --p-custom-axes "$k" \
    --o-visualization core-metrics-results/unweighted_unifrac_emperor-"$k".qzv \

  qiime emperor plot \
  --i-pcoa core-metrics-results/bray_curtis_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results/bray-curtis-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results/jaccard_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results/jaccard-"$k".qzv

  qiime emperor plot \
  --i-pcoa core-metrics-results/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization core-metrics-results/weighted_unifrac-"$k".qzv

  ) &

done;
wait


# Alpha rarefaction plotting
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-min-depth $MIN_DEPTH \
  --p-max-depth $SAMPLING_DEPTH \
  --p-steps 200 \
  --p-iterations 10 \
  --m-metadata-file "$METADATA" \
  --o-visualization alpha-rarefaction.qzv

conda deactivate

#### filtered analysis (remove groups of 1 and 2 samples) ####
##### Define Variables #####

METADATA="/data/gencore/analysis_projects/completed_projects/5920921_Sandrin_16S/metadata-filtered.txt"
DIRECTORY="/data/gencore/analysis_projects/completed_projects/5920921_Sandrin_16S/qiime2"
SAMPLING_DEPTH=10759 #typically, lowest feature count from samples, found in table.qzv
MIN_DEPTH=100 #something lower, to get a good curve from min_depth to sampling_depth

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
#classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2022.2/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2022.2/gg-13-8-99-515-806-nb-classifier.qza"

#put all categorical columns here for unifrac plots by category
categories="DiseaseState
            Sex
            Race
            Diagnosis"

numericals="Age"

##### Run Analysis #####

source activate /data/biocore/programs/conda-envs/qiime2-2022.2

cd "$DIRECTORY"

qiime feature-table filter-samples \
  --i-table table.qza \
  --m-metadata-file $METADATA \
  --o-filtered-table index-filtered-table.qza

# Taxonomic analysis with Silva, if selected
for i in ${!classifiers[@]};
do

(echo "$i"

qiime feature-classifier classify-sklearn \
  --i-classifier ${classifiers["$i"]} \
  --i-reads rep-seqs.qza \
  --o-classification "$i"-filtered-taxonomy.qza

qiime metadata tabulate \
  --m-input-file "$i"-filtered-taxonomy.qza \
  --o-visualization "$i"-filtered-taxonomy.qzv

qiime taxa barplot \
  --i-table index-filtered-table.qza \
  --i-taxonomy "$i"-filtered-taxonomy.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization "$i"-filtered-taxa-bar-plots.qzv

qiime taxa collapse \
    --i-table filtered-table.qza \
    --i-taxonomy "$i"-filtered-taxonomy.qza \
  --o-collapsed-table "$i"-filtered-level7_table.qza \
    --p-level 7

qiime feature-table relative-frequency \
--i-table "$i"-filtered-level7_table.qza \
--o-relative-frequency-table "$i"-filtered-rel-level7_table.qza

mkdir "$i"-filtered-rel-table7
qiime tools export \
--input-path "$i"-filtered-rel-level7_table.qza \
--output-path "$i"-filtered-rel-table7
#this step will probably produce a python error but it succeeds
cd "$i"-filtered-rel-table7
biom convert -i filtered-feature-table.biom -o "$i"-filtered-rel-level7-table.tsv --to-tsv
cd ..) &

done;
wait

### ANCOMBC at class level
module load mamba/latest

DIRECTORY="/data/gencore/analysis_projects/6209010_Keaton_Meta/5920921_Sandrin_16S"
FASTQ_DIR="$DIRECTORY"/fastq
OUTPUT_DIR="$DIRECTORY"/qiime2

METADATA="$DIRECTORY"/metadata-filtered.txt

source activate /data/biocore/programs/conda-envs/qiime2-2022.2

qiime feature-classifier classify-sklearn \
  --i-classifier "/data/biocore/qiime2_classifiers/qiime-2022.2/gg-13-8-99-515-806-nb-classifier.qza" \
  --i-reads rep-seqs.qza \
  --o-classification gg-filtered-taxonomy.qza

qiime taxa collapse \
  --i-table index-filtered-table.qza \
  --i-taxonomy gg-filtered-taxonomy.qza \
  --o-collapsed-table gg-filtered-level3_table.qza \
  --p-level 3

qiime taxa collapse \
  --i-table index-filtered-table.qza \
  --i-taxonomy gg-filtered-taxonomy.qza \
  --o-collapsed-table gg-filtered-level5_table.qza \
  --p-level 5

qiime taxa collapse \
  --i-table index-filtered-table.qza \
  --i-taxonomy gg-filtered-taxonomy.qza \
  --o-collapsed-table gg-filtered-level7_table.qza \
  --p-level 7

source deactivate
source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2024.10-try2

qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/gg-filtered-level3_table.qza \
  --m-metadata-file "$METADATA" \
  --p-formula DiagnosisType \
  --p-reference-levels DiagnosisType::Healthy \
  --o-differentials "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level3.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level3.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-diagnosistype.healthy-gg-filtered-level3-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level3.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-diagnosistype.healthy-gg-filtered-level3.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level3.qza \
  --output-path "$OUTPUT_DIR"/ancombc-filtered-level3

qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/gg-filtered-level5_table.qza \
  --m-metadata-file "$METADATA" \
  --p-formula DiagnosisType \
  --p-reference-levels DiagnosisType::Healthy \
  --o-differentials "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level5.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level5.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-diagnosistype.healthy-gg-filtered-level5-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level5.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-diagnosistype.healthy-gg-filtered-level5.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level5.qza \
  --output-path "$OUTPUT_DIR"/ancombc-filtered-level5

m
qiime composition ancombc \
  --i-table "$OUTPUT_DIR"/gg-filtered-level7_table.qza \
  --m-metadata-file "$METADATA" \
  --p-formula DiagnosisType \
  --p-reference-levels DiagnosisType::Healthy \
  --o-differentials "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level7.qza

qiime composition da-barplot \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level7.qza \
  --p-effect-size-label 'lfc' \
  --p-feature-id-label 'id' \
  --p-error-label 'se' \
  --p-significance-label 'q_val' \
  --p-significance-threshold 0.05 \
  --p-effect-size-threshold 0.0 \
  --p-level-delimiter ';' \
  --o-visualization "$OUTPUT_DIR"/dabarplot-ancombc-diagnosistype.healthy-gg-filtered-level7-qval.05.qza

qiime composition tabulate \
  --i-data "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level7.qza \
  --o-visualization "$OUTPUT_DIR"/table-ancombc-diagnosistype.healthy-gg-filtered-level7.qzv

qiime tools export \
  --input-path "$OUTPUT_DIR"/ancombc-diagnosistype.healthy-gg-filtered-level7.qza \
  --output-path "$OUTPUT_DIR"/ancombc-filtered-level7
