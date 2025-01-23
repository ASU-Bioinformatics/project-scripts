#!/bin/bash

##### trying to get abundance tables #####

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

cd /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/coassembly-contig-annotations/htseqcounts

# change the sample prefix to any sample in the dataset
awk -F\\ '{print $1}' P5_SQP_full.htseqcounts.txt > transcripts.txt
sed -i "1s/^/transcript_id \n/" transcripts.txt

for i in $(find ./ -type f -name "*htseqcounts.txt" | while read F; do basename $F; done | cut -d "." -f 1 )
do
  echo "$i"
  awk -F 't' '{print $NF}' "$i.htseqcounts.txt" > "$i.counts-only.txt"
  sed -i "1s/^/$i \n/" "$i.counts-only.txt"
done

paste transcripts.txt $(find ./ -type f -name "*.counts-only.txt" | sort | uniq ) > gene.count.matrix.tsv
