#!/bin/bash
#SBATCH -c 1
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-1:00
#SBATCH --mem=4G
#SBATCH -o slurm.%A.bedtools-cov.out
#SBATCH -e slurm.%A.bedtools-cov.err

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/pureclip-env

mdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/mouse"
vdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/vaccinia"

#sort -k1,1 -k2,2n \
#     "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

#sort -k1,1 -k2,2n \
#     "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

#sort -k1,1 -k2,2n \
#     "$mdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$mdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

# these parameters create a histogram showing features in the genome that are covered by 90% of the pureClip peak
# since features can be fairly long, it may be only a small percent of the genomic feature itself
# for the contrast between IP and Input vaccinia, there are 210 features bound by the protein above noise level
# I could find out the more precise probabilities from the original bed file.
#bedtools coverage \
#  -b "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
#  -a /data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.gff \
#  -hist -F 0.90 -s -sorted -header \
#  > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.txt

#cat 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.txt \
#| tr ' ' '_' | awk '$10 ~ /1/ { print $0 }' \
#> 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.present.txt

bedtools coverage \
  -b "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.114.sorted.gff3 \
  -hist -s -sorted -header \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist.txt

cat "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist.txt \
| tr ' ' '_' | awk '$10 ~ /1/ { print $0 }' \
> "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist.present.txt

bedtools coverage \
  -b "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.114.sorted.gff3 \
  -hist -s -F 0.90 -sorted -header \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.txt

cat "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.txt \
| tr ' ' '_' | awk '$10 ~ /1/ { print $0 }' \
> "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.covHist-0.9.present.txt
