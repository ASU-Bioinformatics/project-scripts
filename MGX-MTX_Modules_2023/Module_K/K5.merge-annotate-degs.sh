#!/bin/bash

##### merge DEG lists and add annotations #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.K6out
#SBATCH -e slurm.%j.K6.err
#SBATCH -t 0-1:00
#SBATCH --mem=32G

#### Merge DEG Lists from Three Tools ####
# identify the comparison header names
diffExprDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression"
comparisonFile=$diffExprDir/comparisons.csv
deseqDir=$diffExprDir/automated-deseq
edgerDir=$diffExprDir/automated-edger
noiseqDir=$diffExprDir/automated-noiseq

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

cd $diffExprDir

python /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K5.merge_de.py \
       $comparisonFile $deseqDir $edgerDir $noiseqDir

# in the annotations file column 2 corresponds to the transcript ID used by the DEG lists
# columns 12-23 are functional and taxonomic predictions
# col 12 is a protein name when available (such as a WP_ id)
# col 13 is a description of the protein when col 12 is present
# col 14 is a KEGG ID for the protein
# col 15 is a description of the protein based of the available KEGG ID
# col 16 is a predicted taxonomic classification for the protein
# col 17 is all molecular function GO terms associated with the protein
# col 18 is all cellular component GO terms associated with the protein
# col 19 is all biological process GO terms associated with the protein
# col 20 is the interproscan IDs associated with the protein
# col 21 is the pfam IDs associated with the protein
# col 22 is a series of numbers that I'm not sure about the origin of
# col 23 is the database used to identify the protein
annotations="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/complete-annotations-nocounts.tsv"

LC_ALL=C
LANG=C

sortedAnnotations=$diffExprDir/sorted-annotations.tsv

sort -k2 $annotations > $sortedAnnotations

for i in $(find ./ -maxdepth 1 -type f -name "merged_*.csv" | while read F; do basename $F; done )
do
  echo "$i"
  sort -k1 -t ',' $i | tr ',' '\t' | dos2unix > tmp.txt
  echo $(head tmp.txt)
  join -o 1.1,1.2,1.3,1.4,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.23 -1 1 -2 2 -t $'\t' tmp.txt "$sortedAnnotations" > annotated_${i%csv}txt
done

mkdir unannotated-merged
mv merged* unannotated-merged

mkdir annotated-merged
mv annotated_merged* annotated-merged

cd annotated-merged
for i in $(find ./ -maxdepth 1 -type f -name "annotated_*" | while read F; do basename $F; done )
do
  echo "$i"
  awk -F $'\t' '{ if (($2!="NA" && $3!="NA") || ($2!="NA" && $4!="NA") || ($3!="NA" && $4!="NA")) print }' "$i" > filtered_$i
done

mkdir ../filtered-annotated-merged
mv filtered_annotated_merged_deg.* ../filtered-annotated-merged/

# do this in DNA paired data to get reference KEGG id counts
cat final.out.CAT.predicted_proteins.faa.ko | sort |uniq -c | sort -nr | awk '{$1=$1};1' | tr ' ' '\t' > final.out.CAT.predicted_proteins.faa.keggcounts.txt
for i in $(find ./ -maxdepth 1 -type f -name "*.keggs.tsv" | while read F; do basename $F; done )
do
  echo "$i"
  cat "$i" | sort | uniq -c | sort -nr | awk '{$1=$1};1' | tr ' ' '\t' > ${i%keggs.tsv}keggcounts.tsv
done

for i in $(find ./ -maxdepth 1 -type f -name "filtered_*" | while read F; do basename $F; done )
do
  echo "$i"
  awk -F $'\t' '{ if (($7!="NA") && ($7!="null")) print $7 }' "$i" | sort | uniq -c | sort -nr > keggcounts.RNA.$i
done

for i in $(find ./ -maxdepth 1 -type f -name "filtered_*" | while read F; do basename $F; done )
do
  echo "$i"
  awk -F $'\t' '{ if (($7!="NA") && ($7!="null")) print $7 }' "$i" > kegg.list.RNA.$i
done

kolist=$(find ./ -maxdepth 1 -type f -name "kegg.list.*" | while read F; do basename $F; done)

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/microbe-annotator-env

python /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline/ko_mapper.py \
  -i $kolist \
  -p ko_map --cluster rows
