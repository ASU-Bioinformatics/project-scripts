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

METADATA="/data/gencore/analysis_projects/Marshall/16S/metadata.txt"
DIRECTORY="/data/gencore/analysis_projects/Marshall/16S/qiime2"
SAMPLING_DEPTH=186351 #lowest feature count from samples, found in table.qzv
MIN_DEPTH=100 #something lower, to get a good curve from min_depth to sampling_depth

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2020.8/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2020.8/gg-13-8-99-515-806-nb-classifier.qza"

#put all categorical columns here for unifrac plots by category
categories="BarcodeSequence
            LinkerPrimerSequence
            Description
            Replicate
            Location
            treatment
            Replicates"

##### Run Analysis #####

source activate /data/biocore/programs/conda-envs/qiime2-2021.4

cd "$DIRECTORY"

# Taxonomic analysis with selected classifiers
for i in ${!classifiers[@]};
do

(echo "$i"

qiime feature-classifier classify-sklearn \
  --i-classifier ${classifiers[sv]} \
  --i-reads rep-seqs.qza \
  --o-classification sv-taxonomy.qza

qiime metadata tabulate \
  --m-input-file sv-taxonomy.qza \
  --o-visualization sv-taxonomy.qzv

qiime taxa barplot \
  --i-table table.qza \
  --i-taxonomy "$i"-taxonomy.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization "$i"-taxa-bar-plots.qzv

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
