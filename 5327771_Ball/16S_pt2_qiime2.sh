#!/bin/bash

##### qiime2 pt2: taxonomy and phylogeny #####

#SBATCH -A kawoodbu
#SBATCH -o slurm.%j.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err               # STDERR (%j = JobId)
#SBATCH --mail-type=ALL               # Send a notification when the job starts, stops, or fails
#SBATCH --mail-user=kawoodbu@asu.edu  # send-to address
#SBATCH -t 1-0:00                    # estimated time needed - silva classifier takes a long time
#SBATCH --mem=32G

module purge

DIRECTORY="/data/gencore/analysis_projects/5327771_Ball_16S"
FASTQ_DIR="$DIRECTORY"/fastq
OUTPUT_DIR="$DIRECTORY"/qiime-results
METADATA="$DIRECTORY"/scripts/"metadata.tsv"

# set sampling depth based on min feature count from table.qzv; min depth is typically between 100-1000
SAMPLING_DEPTH=21131
MIN_DEPTH=100

source activate /data/biocore/programs/conda-envs/qiime2-2022.2

#comment out a line to exclude a classifier, or add a new line to include a new classifier
declare -A classifiers
classifiers[sv]="/data/biocore/qiime2_classifiers/qiime-2022.2/silva-138-99-515-806-nb-classifier.qza"
classifiers[gg]="/data/biocore/qiime2_classifiers/qiime-2022.2/gg-13-8-99-515-806-nb-classifier.qza"

categories="SuccessionalStage
            PlantType"

numericals="Lat
            Long"

cd "$OUTPUT_DIR"

# Taxonomic analysis with selected classifiers
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

#standard out
qiime diversity core-metrics-phylogenetic \
  --i-phylogeny rooted-tree.qza \
  --i-table table.qza \
  --p-sampling-depth $SAMPLING_DEPTH \
  --m-metadata-file "$METADATA" \
  --output-dir "core-metrics-results-$SAMPLING_DEPTH"

qiime diversity alpha-group-significance \
  --i-alpha-diversity "core-metrics-results-$SAMPLING_DEPTH"/faith_pd_vector.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/faith-pd-group-significance.qzv

qiime diversity alpha-group-significance \
  --i-alpha-diversity "core-metrics-results-$SAMPLING_DEPTH"/evenness_vector.qza \
  --m-metadata-file "$METADATA" \
  --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/evenness-group-significance.qzv

#beta group significance, for each categorical column, standard and decon outs
for j in $categories;
do

  (echo "$j"

  qiime diversity beta-group-significance \
    --i-distance-matrix "core-metrics-results-$SAMPLING_DEPTH"/unweighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/unweighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix "core-metrics-results-$SAMPLING_DEPTH"/weighted_unifrac_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/weighted_unifrac_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix "core-metrics-results-$SAMPLING_DEPTH"/bray_curtis_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/bray_curtis_"$j"_significance.qzv \
    --p-pairwise

  qiime diversity beta-group-significance \
    --i-distance-matrix "core-metrics-results-$SAMPLING_DEPTH"/jaccard_distance_matrix.qza \
    --m-metadata-file "$METADATA" \
    --m-metadata-column "$j" \
    --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/jaccard_"$j"_significance.qzv \
    --p-pairwise

  ) &

done;
wait

for k in $numericals;
do

  echo "$k"

  qiime emperor plot \
    --i-pcoa "core-metrics-results-$SAMPLING_DEPTH"/unweighted_unifrac_pcoa_results.qza \
    --m-metadata-file "$METADATA" \
    --p-custom-axes "$k" \
    --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/unweighted_unifrac_emperor-"$k".qzv

  qiime emperor plot \
  --i-pcoa "core-metrics-results-$SAMPLING_DEPTH"/bray_curtis_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/bray-curtis-"$k".qzv

  qiime emperor plot \
  --i-pcoa "core-metrics-results-$SAMPLING_DEPTH"/jaccard_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/jaccard-"$k".qzv

  qiime emperor plot \
  --i-pcoa "core-metrics-results-$SAMPLING_DEPTH"/weighted_unifrac_pcoa_results.qza \
  --m-metadata-file "$METADATA" \
  --p-custom-axes "$k" \
  --o-visualization "core-metrics-results-$SAMPLING_DEPTH"/weighted_unifrac-"$k".qzv


done;
wait

# Alpha rarefaction plotting standard out
qiime diversity alpha-rarefaction \
  --i-table table.qza \
  --i-phylogeny rooted-tree.qza \
  --p-min-depth $MIN_DEPTH \
  --p-max-depth $SAMPLING_DEPTH \
  --p-steps 200 \
  --p-iterations 10 \
  --m-metadata-file "$METADATA" \
  --o-visualization "alpha-rarefaction.$SAMPLING_DEPTH.qzv"

#leave qiime2 environment so that python3.9 can be used to generate sunburst plots
#source deactivate

#for i in ${!classifiers[@]};
#do
#  (python /data/gencore/shared_scripts/sunburst_csvs.py "$i" "$STD_DIR"
#  python /data/gencore/shared_scripts/sunburst_csvs.py "$i" "$DECON_DIR") &
#done

#wait

#for i in ${!classifiers[@]};
#do

#(echo "$i"
#qiime taxa collapse \
#    --i-table table.qza \
#    --i-taxonomy "$OUTPUT_DIR"/"$i"-taxonomy.qza \
#    --o-collapsed-table "$OUTPUT_DIR"/"$i"-level5_table.qza \
#    --p-level 5

#qiime feature-table relative-frequency \
#--i-table "$OUTPUT_DIR"/"$i"-level5_table.qza \
#--o-relative-frequency-table "$OUTPUT_DIR"/"$i"-rel-level5_table.qza

#mkdir "$OUTPUT_DIR"/"$i"-rel-table5
#qiime tools export \
#--input-path "$OUTPUT_DIR"/"$i"-rel-level5_table.qza \
#--output-path "$OUTPUT_DIR"/"$i"-rel-table5
#this step will probably produce a python error but it succeeds
#cd "$OUTPUT_DIR"/"$i"-rel-table5
#biom convert -i feature-table.biom -o "$OUTPUT_DIR"/"$i"-rel-level5-table.tsv --to-tsv
#cd ..) &

#done;
