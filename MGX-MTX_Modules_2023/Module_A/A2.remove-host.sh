#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 2-0:00
#SBATCH -c 4
#SBATCH --mem=24G

sid="$1"
assembly="$2"
fqDir="$3"
outDir="$4"

#sid="MB-055"
#assembly="/data/gencore/databases/reference_genomes/human/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome"
#fqDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/fastq/adapter-trimmed/fastq/"
#outDir="/data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq"

echo "bowtie alignment on $assembly for sample $sid"

module purge
module load mamba/latest

source activate /data/biocore/programs/conda-envs/bowtie2-env/
export PERL5LIB=/data/biocore/programs/conda-envs/bowtie2-env/bin/perl

cd "$outDir"

# use SME for untrimmed data, SQP for trimmed data
# also check lane ID. It's usually L001.

bowtie2 --verbose -p 16 -x "$assembly" \
        -1 "$fqDir"/"$sid"_SQP_L001_R1_001.fastq.gz \
        -2 "$fqDir"/"$sid"_SQP_L001_R2_001.fastq.gz \
        --un-conc-gz "$outDir"/"$sid"_nohost_SQP_L001_R%_001.fastq.gz \
        > "$outDir"/"$sid"_SQP_host_mapped_and_unmapped.sam

#echo "full.sam file is complete"

tail -n+5 "$sid"_SQP_host_mapped_and_unmapped.sam > "$sid"_SQP_full.headed.host_mapped_and_unmapped.sam

echo "full.headed.sam file is complete"

# obtain coverage stats and RPKM values for all alignments

source deactivate

module load bbmap-39.01-gcc-12.1.0

pileup.sh in="$outDir"/"$sid"_SQP_full.headed.host_mapped_and_unmapped.sam \
          out="$outDir"/"$sid"_covstats.txt \
          rpkm="$outDir"/"$sid"_rpkm.txt


module load fastqc-0.11.9-gcc-12.1.0

fastqc -t 4 "$outDir"/"$sid"_nohost_SQP_L001_R*_001.fastq.gz
mv "$outDir"/"$sid"_nohost_SQP_L001_R*_001_fastqc* "$outDir"/qc

mv "$outDir"/"$sid"_*.sam "$outDir"/alignments
mv "$outDir"/"$sid"_*.txt "$outDir"/alignments
