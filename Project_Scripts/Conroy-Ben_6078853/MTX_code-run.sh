# Module K, script K1, subscript DESeq2 analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script deseq2 \
  --directory /data/gencore/analysis_projects/6078853_Otak_RNA/differential-expression \
  --geneMatrix /data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/htseqcounts/transcript.count.matrix.tsv \
  --comparisons /data/gencore/analysis_projects/6078853_Otak_RNA/comparisons.csv

# Module K, script K1, subscript edgeR analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script edger \
  --directory /data/gencore/analysis_projects/6078853_Otak_RNA/differential-expression \
  --geneMatrix /data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/htseqcounts/transcript.count.matrix.tsv \
  --comparisons /data/gencore/analysis_projects/6078853_Otak_RNA/comparisons.csv

# Module K, script K1, subscript NOIseq analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script NOIseq \
  --directory /data/gencore/analysis_projects/6078853_Otak_RNA/differential-expression \
  --geneMatrix /data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/htseqcounts/transcript.count.matrix.tsv \
  --comparisons /data/gencore/analysis_projects/6078853_Otak_RNA/comparisons.csv

# make annotation file
gff="/data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/out.CAT.predicted_proteins.gff"
newgff="nocomments.out.CAT.predicted_proteins.gff"
sortedgff="sorted.nocomments.out.CAT.predicted_proteins.gff"
annotations="suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"

anns="/data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/all.out.CAT.predicted_proteins.faa.annot"
sortedanns="sorted.all.out.CAT.predicted_proteins.faa.annot"

matrix="/data/gencore/analysis_projects/6078853_Otak_RNA/coassembly-contig-annotations/htseqcounts/transcript.count.matrix.tsv"
tabbedGFF="transcript.count.matrix_tabbedgff.tsv"
sortedTabbedGFF="sorted.transcript.count.matrix_tabbedgff.tsv"

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

sort -k 8 "$tabbedGFF" > "$sortedTabbedGFF"
join -o auto -a1 -e "null" -1 8 -2 1 -t $'\t' "$sortedTabbedGFF" "$sortedanns" > annotated.transcript.count.matrix.tsv
# the value for sorting the tabbed gff file should be two more than the number of samples in the matrix

sort -k 2 "$annotations" | tr ' ' '\t' > "$sortedAnnotations"
join -o auto -a1 -e "null" -1 2 -2 1 -t $'\t' "$sortedAnnotations" "$sortedanns" > $geneInfo

# the geneInfo file is good for annotating paired MTX data!

# add header line to geneInfo file
sed -i '1i\ORFid\tGFFid\tContigID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tgffAttributes\tproteinID\tproteinDesc\tkeggID\tkeggDesc\ttaxonomy\tGO:MF\tGO:CC\tGO:BP\tIPRid\tpfamID\tECnum\tdatabase' $geneInfo

# annotate differential files
annotations="complete-annotations-nocounts.tsv"

LC_ALL=C
LANG=C

sortedAnnotations="sorted-annotations.tsv"

sort -k2 $annotations > $sortedAnnotations

for i in $(find ./ -maxdepth 1 -type f -name "merged_*.csv" | while read F; do basename $F; done )
do
  echo "$i"
  sort -k1 -t ',' $i | tr ',' '\t' | dos2unix > tmp.txt
  echo $(head tmp.txt)
  join -o 1.1,1.2,1.3,1.4,2.12,2.13,2.14,2.15,2.16,2.17,2.18,2.19,2.20,2.21,2.22,2.23 -1 1 -2 2 -t $'\t' tmp.txt "$sortedAnnotations" > annotated_${i%csv}txt
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

# sort and merge NCBI viral covstats files

{ head -n 1; sort; } < 49_covstats.txt | awk -F "\t" '{ print $1 }' > viral-refseq-avgfoldcov.txt
for i in $(find ./ -maxdepth 1 -type f -name "*covstats.txt" | while read F; do basename $F; done )
do
  echo "${i%_covstats.txt}"
  { head -n 1; sort; } < "$i" | awk -F "\t" '{ print $2 }' > "${i%_covstats.txt}".tmp.txt
  sed -i "1s/.*/${i%_covstats.txt}/" "${i%_covstats.txt}".tmp.txt
  paste viral-refseq-avgfoldcov.txt "${i%_covstats.txt}".tmp.txt > tmp.txt
  cat tmp.txt > viral-refseq-avgfoldcov.txt
done

# stats for alignment to coassembly - run for all 6 samples
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/DX.alignment-stats-1sample.sh 49

# determine kegg/go terms from the annotated count matrix
# the column numbers are based on the number of samples, since this is the count matrix. This is for a run with 6 samples.
awk -F "\t" '{ print $20 }' annotated.transcript.count.matrix.tsv | sort | uniq > uniq.keggs.txt
wc -l uniq.keggs.txt #6,696, -2 for NA and null makes 6,694

awk -F "\t" '{ print $23 }' annotated.transcript.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goMF.txt
wc -l uniq.goMF.txt #2252, -2 for NA and null makes 2,250

awk -F "\t" '{ print $24 }' annotated.transcript.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goCC.txt
wc -l uniq.goCC.txt #713, -2 for NA and null makes 711

awk -F "\t" '{ print $25 }' annotated.transcript.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goBP.txt
wc -l uniq.goBP.txt #2971, -2 for NA and null makes 2,969

# night of Jan 18:
# need to incorporate flagstats and make image, after those runs complete
# create GFF file for the RefSeq viral alignments


# rerun bowtie alignment for sample 60, there appears to be some error in the sam/bam files
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D3.bowtie-align.wrapper.sh \
  --fastqDir /data/gencore/sftp/otakuye_conroyben/6078853_Metatranscriptomics/2.trimmed-filtered_reads/fastq \
  --assembly /data/gencore/sftp/otakuye_conroyben/6078853_Metatranscriptomics/4.metatranscriptome-coassembly/bowtie-index/coassembly \
  --outDir /data/gencore/analysis_projects/6078853_Otak_RNA/megahit_alignments/60_redo \
  --list "60"

# once script D3 is run, all the 60 files can be added to the sftp folder and the excel file can be updated with the 60 data
# then alignment stats image can be generated, and added to the overview document
# and then both MGX and MTX can be returned!
