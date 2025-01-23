#!/bin/bash

##### humann3 test #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.C2.out
#SBATCH -o slurm.%j.C2.err
#SBATCH -t 0-00:30
#SBATCH -c 1
#SBATCH --mem=8G

#### I usually do this one manually to verify the steps are working - it's really fast anyway

module load mamba/latest
source activate /data/biocore/programs/conda-envs/humann-env/

hmDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann90-paired_pipeline"
cd "$hmDir"

for i in $(find ./ -maxdepth 1 -type f -name "*genefamilies.tsv" | while read F; do basename $F; done)
do
  humann_renorm_table --input "$hmDir"/"$i" \
                    --output "$hmDir"/${i%.tsv}"-cpm.tsv" \
                    --units cpm --update-snames
done;

for i in $(find ./ -maxdepth 1 -type f -name "*pathabundance.tsv" | while read F; do basename $F; done)
do
  humann_renorm_table --input "$hmDir"/"$i" \
                    --output "$hmDir"/${i%.tsv}"-copm.tsv" \
                    --units cpm --update-snames
done;

humann_join_tables -i "$hmDir" -o normalized_pathabundance.tsv --file_name pathabundance-copm.tsv

humann_join_tables -i "$hmDir" -o normalized_genefamilies.tsv --file_name genefamilies-cpm.tsv

# create unstratified versions of normalized joined tables, for reducing input to MaAsLin
humann_split_stratified_table --input normalized_pathabundance.tsv --output "$hmDir"

humann_split_stratified_table --input normalized_genefamilies.tsv --output "$hmDir"

# merge metaphlan abundance profiles - not relevant for the RNA part of paired transcriptomics

merge_metaphlan_tables.py ./*/*_metaphlan_bugs_list.tsv > merged_metaphlan_table.txt

grep -E "s__|clade" merged_metaphlan_table.txt \
| grep -v "t__" \
| sed "s/^.*|//g" \
| sed "s/SRS[0-9]*-//g" \
> merged_metaphlan_table_species.txt

grep -E "g__|clade" merged_metaphlan_table.txt \
| grep -v "s__" \
| sed "s/^.*|//g" \
| sed "s/SRS[0-9]*-//g" \
> merged_metaphlan_table_genus.txt

grep -E "f__|clade" merged_metaphlan_table.txt \
| grep -v "g__" \
| sed "s/^.*|//g" \
| sed "s/SRS[0-9]*-//g" \
> merged_metaphlan_table_family.txt

grep -E "o__|clade" merged_metaphlan_table.txt \
| grep -v "f__" \
| sed "s/^.*|//g" \
| sed "s/SRS[0-9]*-//g" \
> merged_metaphlan_table_order.txt

# this typically doesn't give good results
hclust2.py \
-i merged_metaphlan_table_species.txt \
-o metaphlan4_abundance_heatmap_species.png \
--f_dist_f braycurtis \
--s_dist_f braycurtis \
--cell_aspect_ratio 0.5 \
--flabel_size 2 --slabel_size 2 \
--max_flabel_len 20 --max_slabel_len 20 \
--dpi 700 --image_size 14

#--log_scale \
#--minv 0.1 \

# generate metaphlan heatmap with metaphlan-heatmap.R

# organize files for data return
# don't do this for paired samples until the paired RNA runs are done!
mkdir -p pathways/pathway-abundance/samples-normalized
mkdir -p pathways/pathway-abundance/samples-raw
mkdir -p pathways/pathway-coverage
mkdir -p pathways/differential-pathways
mkdir -p genefamilies/samples-normalized
mkdir -p genefamilies/samples-raw
mkdir -p humann-temp
mkdir -p taxonomic-classifications/samples-metaphlan-alignments
mkdir -p taxonomic-classifications/samples-taxonomic-data

mv *genefamilies.tsv genefamilies/samples-raw
mv *genefamilies-cpm.tsv genefamilies/samples-normalized
mv normalized_genefamilies* genefamilies

mv *pathcoverage.tsv pathways/pathway-coverage
mv normalized_pathabundance* pathways/pathway-abundance
mv *pathabundance.tsv pathways/pathway-abundance/samples-raw
mv *pathabundance-copm.tsv pathways/pathway-abundance/samples-normalized

mv *_humann_temp/*bugs* taxonomic-classifications/samples-taxonomic-data
mv *_humann_temp/*bowtie2.txt taxonomic-classifications/samples-metaphlan-alignments
mv *metaphlan* taxonomic-classifications

mv *_humann_temp humann-temp
