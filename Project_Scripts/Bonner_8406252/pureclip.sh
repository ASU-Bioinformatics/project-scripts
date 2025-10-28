#!/bin/bash
#SBATCH -c 6
#SBATCH -p general
#SBATCH -q public
#SBATCH -t 7-0:00
#SBATCH --mem=128G
#SBATCH -o slurm.%A.pureclip.out
#SBATCH -e slurm.%A.pureclip.err

# and I'm not sure how much time and memory this will take
# I'm going to try the first mouse with 7 days and 256G and see how that goes;
# simultaneously I'll try the first vaccinia with 2 days and 128G
# retrying vaccinia with 7 days and 256G, but increasing thread parameter to 128

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/pureclip-env

mdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/mouse"
vdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/vaccinia"

# code for mouse

echo "Analyzing Mouse Mock Subsampled Small"
pureclip \
  -i "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -bai "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai \
  -g /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  -nt 256 \
  -o "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.runNC.bed \
  -or "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.runNC.merged.bed \
  -ibam "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -ibai "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai

#pureclip \
#  -i "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam \
#  -bai "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai \
#  -g /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa \
#  -nt 256 \
#  -o "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.run2.bed \
#  -or "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#  -ibam "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam \
#  -ibai "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai

# code for vaccinia

#pureclip \
#  -i "$vdir"/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
#  -bai "$vdir"/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai \
#  -g /data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.fa \
#  -nt 128 \
#  -o "$vdir"/Mock.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.bed \
#  -ibam "$vdir"/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
#  -ibai "$vdir"/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai

#pureclip \
#  -i "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
#  -bai "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai \
#  -g /data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.fa \
#  -nt 256 \
#  -o "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.bed \
#  -or "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#  -ibam "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
#  -ibai "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai
