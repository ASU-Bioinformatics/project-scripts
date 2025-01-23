## Analysis for Variant Calling of Two E. coli Genomes

### Alignment to Reference

Since I'm repeating earlier analysis types and comparing these both to the reference genome, I'm starting by aligning the reads representing each genome to the specified *E. coli* genome, HG738867.1.

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/A1.alignment.sh \
              -a /data/gencore/analysis_projects/7484728_Rajeev/align \
              -r /data/gencore/databases/reference_genomes/ecoli \
              -g HG738867.1.genomic.fasta \
              -f /data/gencore/analysis_projects/7484728_Rajeev/fastq
```
This alignment looked exceptionally good! About 97% properly paired mapping for both samples.

### Variant Calling
These three tools focus on different sizes of variants, so they are all worth running (except when genotyping; in that case, only gatk-gvcf should be run).

#### GATK HaplotypeCaller
I'm using the default filtering conditions from this script!

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B1.gatk.sh \
        --inputDir /data/gencore/analysis_projects/7484728_Rajeev/align \
        --ref /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        --outDir /data/gencore/analysis_projects/7484728_Rajeev/variant-calls
```

#### GRIDSS2

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B2.gridss2.sh \
        --inputDir /data/gencore/analysis_projects/7484728_Rajeev/align \
        --ref /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        --outDir /data/gencore/analysis_projects/7484728_Rajeev/variant-calls \
        --scriptDir /data/biocore/programs/gridss-2.13.2/
```

Looking at the GRIDSS2 documentation, it seems like joint calling might be a better option, so I am trying that also. There seemed to me to be a lot of potential rearrangment events in 3361 when I looked at it, and having more certainty about whether those are shared variants would be helpful for comparing the two lab strains directly. On the other hand, the joint calling is throwing errors while the single sample calling wasn't, and there isn't good support for deciphering that on their GitHub.

```
scriptDir=/data/biocore/programs/gridss-2.13.2
outDir=/data/gencore/analysis_projects/7484728_Rajeev/variant-calls
inputDir=/data/gencore/analysis_projects/7484728_Rajeev/align
ref=/data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta

$scriptDir/scripts/gridss \
  -t 8 -r $ref \
  -j $scriptDir/gridss-2.13.2-gridss-jar-with-dependencies.jar \
  -o $outDir/joint.gridss.vcf \
  $inputDir/3361.sorted.rdgrp.mrkd.bam $inputDir/3361.sorted.rdgrp.mrkd.bam

vcftools --vcf $outDir/joint.gridss.vcf \
         --out $outDir/joint.gridss.hardfiltered.vcf \
         --remove-filtered-all --recode

rename 'vcf.recode.vcf' 'vcf' *
```

#### Delly

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B3.delly.sh \
        --inputDir /data/gencore/analysis_projects/7484728_Rajeev/align \
        --ref /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        --outDir /data/gencore/analysis_projects/7484728_Rajeev/variant-calls
```

#### Merge VCF Files

Since there are three VCF files for each sample (one from each tool), we can use bcftools concat to merge them into a single report output file. The GRIDSS2 output files do have a much different goal and reporting style than the other two, so it may be preferable just to merge Delly and GATK and report GRIDSS2 separately.

```
module load bcftools-1.14-gcc-11.2.0
bcftools concat 3361.filtered.vcf 3361.gridss.hardfiltered.vcf 3361.delly.hardfiltered.vcf | bcftools sort > 3361.merged.vcf
bcftools concat 3361.filtered.vcf 3361.delly.hardfiltered.vcf | bcftools sort > 3361.gatk.delly.merged.vcf

bcftools concat 3362.filtered.vcf 3362.gridss.hardfiltered.vcf 3362.delly.hardfiltered.vcf | bcftools sort > 3362.merged.vcf
bcftools concat 3362.filtered.vcf 3362.delly.hardfiltered.vcf | bcftools sort > 3362.gatk.delly.merged.vcf
```

## So Actually...

Rajeev wanted analysis similar to what we ran before, so I'm going to just copy the code I used last time for this set of samples. Sigh.

The alignment results from the first try are still perfect at least! I just need to run a gvcf variant call with a haploid perspective.

For the variant calling:
```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B4.gatk-gvcf.sh \
        -i /data/gencore/analysis_projects/7484728_Rajeev/align \
        -r /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        -o /data/gencore/analysis_projects/7484728_Rajeev/align/variant-calls-gvcf-haploid \
        -n 'misra.all'

mv /data/gencore/analysis_projects/7484728_Rajeev/align/variant-calls-gvcf-haploid /data/gencore/analysis_projects/7484728_Rajeev/variant-calls-gvcf-haploid
```

For the SnpEff annotation:

```
sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/7484728_Rajeev/variant-calls-gvcf-haploid \
  -r Escherichia_coli_str_k_12_substr_mc4100 \
  -o /data/gencore/analysis_projects/7484728_Rajeev/variant-calls-gvcf-haploid \
  -s /data/biocore/programs/snpEff
```

For population comparisons after annotation, the R file included in this folder provides a solid workflow, but relatedness isn't really applicable with only two samples.

For the data return, I changed the prefix 'misra.all' to 'both' as it seemed more accurate.

Since Rajeev specifically wanted to see differences between the two species, I'm filtering the VCF file to retain only locations where the two genotypes are different (one matching the reference, one a SNP)

```
bcftools view -e 'GT[0]="R" && GT[1]="R" || GT[0]="A" && GT[1]="A"' both.genotyped.gatk.g.snpeff.vcf > differences.genotyped.gatk.g.snpeff.vcf
```
