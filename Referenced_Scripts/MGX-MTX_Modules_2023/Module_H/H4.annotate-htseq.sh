#!/bin/bash

##### htseq annotation test #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.H4.out
#SBATCH -e slurm.%j.H4.err
#SBATCH -t 0-1:00
#SBATCH --mem=32G

gff="../contig-redo-cat/out.CAT.predicted_proteins.gff"
newgff="nocomments.out.CAT.predicted_proteins.gff"
sortedgff="sorted.nocomments.out.CAT.predicted_proteins.gff"
annotations="suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"

anns="all.out.CAT.predicted_proteins.faa.annot"
sortedanns="sorted.all.out.CAT.predicted_proteins.faa.annot"

matrix="/data/gencore/sftp/otakuye_conroyben/6078853_Metagenomics/8.abundance-profiling/gene.count.matrix.tsv"
tabbedGFF="gene.count.matrix_tabbedgff.tsv"
sortedTabbedGFF="sorted.gene.count.matrix_tabbedgff.tsv"

sortedAnnotations="sorted.suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"
geneInfo="complete-annotations-nocounts.tsv"

sed '/^#/ d' < "$gff" > "$newgff"

# sort the files
LC_ALL=C
LANG=C
sort -k 9 "$newgff" > "$sortedgff"
sort -k 1 "$anns" > "$sortedanns"

awk '{ split($9, a, ";"); print a[1], $1, $2, $3, $4, $5, $6, $7, $8, $9; }' "$sortedgff" \
| awk '{ split($1, a, "="); print a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10; }' \
| awk '{ split($1, a, "_"); print $1, $2 "_" a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10; }' > "$annotations"

LC_ALL=C
LANG=C
join -o auto -a1 -e "null" <(sort "$matrix") <(sort "$annotations") | tr ' ' '\t' > "$tabbedGFF"

LC_ALL=C
LANG=C

sort -k 16 "$tabbedGFF" > "$sortedTabbedGFF"
join -o auto -a1 -e "null" -1 16 -2 1 -t $'\t' "$sortedTabbedGFF" "$sortedanns" > annotated.gene.count.matrix.tsv
# the value for sorting the tabbed gff file should be two more than the number of samples in the matrix

sort -k 2 "$annotations" | tr ' ' '\t' > "$sortedAnnotations"
join -o auto -a1 -e "null" -1 2 -2 1 -t $'\t' "$sortedAnnotations" "$sortedanns" > $geneInfo

# the geneInfo file is good for annotating paired MTX data!

# add header line to geneInfo file
sed -i '1i\ORFid\tGFFid\tContigID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tgffAttributes\tproteinID\tproteinDesc\tkeggID\tkeggDesc\ttaxonomy\tGO:MF\tGO:CC\tGO:BP\tIPRid\tpfamID\tECnum\tdatabase' $geneInfo

# after annotating, use R to create a normalized matrix
# R code: df[-1] <- sapply(df[-1], prop.table) * 100

#### this is an inefficient method! better than grep, which in turn is better than python, but still very slow ####
#### do not use! but keep as a reminder of poor programming lol ####
# even with awk this is a slow script, split gene matrix into 20 chunks to parallelize the annotation
#split gene.count.matrix.tsv -n 20 gene_counts_split.tsv_ -d

#sid=00
#input="gene_counts_split.tsv"_"$sid"

# script for when chromosome info is not included in HTSeq count files
# including chromosome info automatically in HTSeq will reduce the time needed for this script
#while IFS= read -r line
#do
  #echo "$line"
#  tid=$(cut -f1  <<< "$line")
#  #echo "tid is" $tid
#  suffix=$(cut -d "_" -f2 <<< "$tid")
#  #echo "suffix is" $suffix
#  gffLine=$(awk -v tid="ID=$tid" '$9 ~ tid { print; exit; }' "$newgff")
#  #echo $gffLine
#  chr=$(cut -f1 <<< "$gffLine")
#  #echo $chr
#  annoCheck="$chr"_"$suffix"
#  #echo $annoCheck
#  annotation=$(awk -v anc="$annoCheck" '$1 ~ anc { print; exit; }' "$anns")
#  #echo $annotation
#  newline="$line"" ""$annotation"
#  echo "$newline"
#done < "$input" > annotated.gene.count.matrix.tsv_"$sid"

# speed test
#tid="1000001_1"
#suffix="1"
#time grep "$tid" "$newgff"
#time awk -v tid="ID=$tid" '$9 ~ tid { print; exit; }' "$newgff"

#awk -v anc="$annoCheck" '$1 ~ anc { print; exit; }' "$anns"
