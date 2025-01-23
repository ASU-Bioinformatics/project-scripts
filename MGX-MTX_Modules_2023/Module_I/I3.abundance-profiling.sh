#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00
#SBATCH -n 4
#SBATCH --mem=64G

module load mamba/latest

source deactivate
source activate /data/biocore/programs/mamba-envs/htseq-env

sid="$1"

echo "$sid htseq quantification using CAT predicted proteins, 8 hours, 1024G memory, 16 threads"

alnDir="$2"

gffFile="$3"

outDir="$4"

htseq-count -f bam -r pos -s no -t CDS -i ID --add-chromosome-info \
            -d "\t" -n 16 --max-reads-in-buffer=30000000000000 \
            "$alnDir"/"$sid".sorted.bam \
            "$gffFile" > "$outDir"/"$sid".htseqcounts.txt

# combine all htseqcounts files to construct raw count matrix

# normalize with deseq2 in R

# take ratio of normalized transcript counts to normalized gene counts
# to get a feel for the relative proportions of transcripts to the functional potential of the population

# compare just the transcript counts with standard DE methods
# but also compare the ratios
