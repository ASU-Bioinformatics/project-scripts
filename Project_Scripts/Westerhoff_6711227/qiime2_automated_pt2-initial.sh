#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 0-8:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=32G

module purge
module load mamba/latest

##### Define Variables #####

METADATA="/data/gencore/analysis_projects/6711227_Partho_16s/metadata.txt"
DIRECTORY="/data/gencore/analysis_projects/6711227_Partho_16s/fastq/qiime2"
SAMPLING_DEPTH=7368 #typically, lowest feature count from samples, found in table.qzv
MIN_DEPTH=100 #something lower, to get a good curve from min_depth to sampling_depth

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2023.7/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2023.7/gg_2022_10_backbone.v4.nb.qza"
#classifiers[18]="/data/biocore/qiime2_classifiers/qiime-2023.7/silva-138-99-nb-classifier.qza"
#classifiers[it]="/data/biocore/qiime2_classifiers/qiime-2022.2/unite-ver9-classifier-99-27-10-2022.qza"

#put all categorical columns here for unifrac plots by category
categories="Soil-type Treatment-GNA-dose Sampling-period"

numericals="Soil-pH Soil-EC Soil-organic-matter Soil-exchangeable-Nitrate Soil-exchangeable-Ammonium Soil-microbial-biomass-C"

##### Run Analysis #####

source activate /data/biocore/programs/mamba-envs/qiime2023.7

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
