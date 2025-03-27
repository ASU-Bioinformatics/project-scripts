## Analysis for Variant Calling of Two E. coli Genomes

### Alignment to Reference

Since I'm repeating earlier analysis types and comparing these both to the reference genome, I'm starting by aligning the reads representing each genome to the specified *E. coli* genome, HG738867.1. (This script must run on Sol.)

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/A1.alignment.sh \
              -a /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align \
              -r /data/gencore/databases/reference_genomes/ecoli \
              -g HG738867.1.genomic.fasta \
              -f /data/gencore/analysis_projects/8087757_Kelly-Rajeev/fastq
```
This alignment looked exceptionally good! About 97% properly paired mapping for both samples.

### Variant Calling

Rajeev wants analysis similar to what we ran before, so I'm going to just copy the code I used last time for this set of samples.

The alignment results are basically perfect! I just need to run a gvcf variant call with a haploid perspective.

For the variant calling:
```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B4.gatk-gvcf.sh \
        -i /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align \
        -r /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        -o /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align/variant-calls-gvcf-haploid \
        -n 'misra.both'

mv /data/gencore/analysis_projects/7484728_Rajeev/align/variant-calls-gvcf-haploid /data/gencore/analysis_projects/7484728_Rajeev/variant-calls-gvcf-haploid
```

### SnpEff Annotation

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8087757_Kelly-Rajeev/variant-calls-gvcf-haploid \
  -r Escherichia_coli_str_k_12_substr_mc4100 \
  -o /data/gencore/analysis_projects/8087757_Kelly-Rajeev/variant-calls-gvcf-haploid \
  -s /data/biocore/programs/snpEff
```

Since Rajeev specifically wanted to see differences between the two species, I'm filtering the VCF file to retain only locations where the two genotypes are different (one matching the reference, one a SNP)

```
bcftools view -e 'GT[0]="R" && GT[1]="R" || GT[0]="A" && GT[1]="A"' both.genotyped.gatk.g.snpeff.vcf > differences.genotyped.gatk.g.snpeff.vcf
```

SnpEff annotation shows the expected SNP in rpoB (double check that it is at base 1295 in each sample), as well as a second point mutation in NQ3-1. There is also the expected additional mutation in rpoC. Identify the base locations of those mutations and send the details to Rajeev as a special note with the results!

### Extract Regions of Interest

For this, I want to create a .bed file containing the chromosome, start position, and end position of the genes of interest, then extract the corresponding fasta from each of the three genomes, and then align them and visualize the alignment in R. I'm particularly interested in regions that appear to be insertions/deletions, seen in the genotyped vcf files; rpoB and rpoC should probably also be highlighted.

HOWEVER this method is DEFINITELY not giving me what I want. It appears to be making a consensus sequence that fills in gaps in the new strain with the bases from the first strain instead of with alignment-type gap indicators, which means that it misses the deleted dnaJ gene in NQ3-1.

```
module load bedtools2-2.30.0-gcc-11.2.0

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/ref-dnaJ.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/ref-dnaJ.fasta

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/ref-dnaJ.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/im3-dnaJ.fasta

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/ref-dnaJ.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/nq3-dnaJ.fasta

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/rpoB.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/ref-rpoB.fasta

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/rpoB.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/im3-rpoB.fasta

bedtools getfasta \
  -fi /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
  -bed /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/rpoB.bed \
  -fo /data/gencore/analysis_projects/8087757_Kelly-Rajeev/beds/nq3-rpoB.fasta
```

In order to do this for the non-reference genomes, I will first need to generate a consensus fasta file from the bam alignment files, which can be done with samtools consensus. I think I will need to specify to include both deletions and insertions, and use the --mark-ins option to preserve coordinate mapping. A first try:

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/D1.samtools-consensus.sh \
  -i /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align/NQ3-1.sorted.rdgrp.mrkd.bam \
  -o /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align/nq3-all.fasta

sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/D1.samtools-consensus.sh \
  -i /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align/IM3.sorted.rdgrp.mrkd.bam \
  -o /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align/im3-all.fasta
```

### Trying Some New Things

I'm not seeing a good way to get the data of interest from the genotyping VCF files, and tbh it doesn't seem to be capturing even the known mutation in rpoB (which is plainly consistent across all reads in the bam files for both samples). So I'm going back to the drawing board and running standard HaplotypeCaller, then of course snpEff annotation afterward.


```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B1.gatk.sh \
        --inputDir /data/gencore/analysis_projects/8087757_Kelly-Rajeev/align \
        --ref /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        --outDir /data/gencore/analysis_projects/8087757_Kelly-Rajeev/variant-calls

sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8087757_Kelly-Rajeev/variant-calls \
  -r Escherichia_coli_str_k_12_substr_mc4100 \
  -o /data/gencore/analysis_projects/8087757_Kelly-Rajeev/variant-calls \
  -s /data/biocore/programs/snpEff
```
