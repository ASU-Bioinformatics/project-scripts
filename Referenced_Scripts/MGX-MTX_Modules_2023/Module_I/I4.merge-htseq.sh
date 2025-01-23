#!/bin/bash

##### get abundance tables #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.htseq-test.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.htseq-test.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00

module load mamba/latest

source deactivate
source activate /data/biocore/programs/mamba-envs/htseq-env

# move all htseqcounts.txt files to a new directory before running
# then cd to that directory

countDir=$1
cd $countDir

# use a random sample to create the transcripts list
rep=$(find ./ -name "*.htseqcounts.txt" | head -n1)

awk -F\\ '{print $1}' $rep > transcripts.txt
sed -i "1s/^/transcript_id \n/" transcripts.txt

for i in $(find ./ -maxdepth 1 -type f -name "*.htseqcounts.txt" | while read F; do basename $F; done | cut -d "." -f 1 )
do
  echo "$i"
  awk -F 't' '{print $NF}' "$i.htseqcounts.txt" > "$i.counts-only.txt"
  sed -i "1s/^/$i \n/" "$i.counts-only.txt"
done

paste transcripts.txt $(find ./ -maxdepth 1 -type f -name "*.counts-only.txt" | sort | uniq ) > transcript.count.matrix.tsv

#rm "*.counts-only.txt"

# not incorporated fully yet! This same strategy can also be used to annotate the deg lists, so maybe I should make it into it's own script...
#geneInfo="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/complete-annotations-nocounts.tsv"
#sortedGeneInfo="/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/sorted-annotations.tsv"

#LC_ALL=C
#LANG=C

#sort -k 2 $geneInfo > $sortedGeneInfo
#sort -k 1 "/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/gene.count.matrix.tsv" > "/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/sorted.gene.count.matrix.tsv"
#join -t $'\t' -a2 -1 2 -2 1 $sortedGeneInfo "/data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/sorted.gene.count.matrix.tsv" > annotated.gene.count.matrix.tsv
