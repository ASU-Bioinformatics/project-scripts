# Variant Calling for Single and Dual Culture Microbes

## Alignment to Reference Genomes

### Single-Culture Alignments

First, I'm aligning the single-culture samples to their appropriate genomes: CAN1102 to *E. coli* and CAN1103 to *Pseudomonas aeruginosa*.

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/A1.alignment.sh \
              -a /data/gencore/analysis_projects/8040442_Yang/align-cft073 \
              -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1 \
              -g CFT073.genomic.fna \
              -l /data/gencore/analysis_projects/8040442_Yang/fastq/CAN1102

sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/A1.alignment.sh \
              -a /data/gencore/analysis_projects/8040442_Yang/align-pao1 \
              -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
              -g GCF_000006765.1_PAO1.genomic.fna \
              -l /data/gencore/analysis_projects/8040442_Yang/fastq/CAN1103
```

### Individual VCF Creation

I don't know how useful this will actually be, since I think the genotyping information will be more valuable for this project, but it could potentially be used as back-up or confirmatory information regarding specific timepoints.

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/B1.gatk.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/align-cft073 \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/CFT073.genomic.fna \
  -o /data/gencore/analysis_projects/8040442_Yang/variants-cft073

sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/B1.gatk.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/align-pao1 \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.genomic.fna \
  -o /data/gencore/analysis_projects/8040442_Yang/variants-pao1
```

### Genotyping VCF Creation

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/B4.gatk-gvcf.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/align-cft073 \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/CFT073.genomic.fna \
  -o /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073 \
  -n 'cft073.genotyping.'

sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/B4.gatk-gvcf.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/align-pao1 \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.genomic.fna \
  -o /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1 \
  -n 'pao1.genotyping.'
```

I'm not sure why, but three of the pao1 samples didn't run and I missed it... here's the code to rerun those:

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/B4.gatk-gvcf.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/align-pao1 \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.genomic.fna \
  -o /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1 \
  -n 'pao1.genotyping.' \
  -l "X29-G117-AG X29-G117-AI X29-G117-AK"
```

### Genotyping VCF Annotation

I'm pretty sure that the gvcf files can't be directly annotated, but the `genotyped` vcfs can because they report SNPs directly rather than regions.

```
java -Xmx4g -jar /data/biocore/programs/snpEff/snpEff.jar \
  -v Escherichia_coli_cft073 /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073/cft073.genotyping.genotyped.gatk.g.vcf \
  > cft073.filtered.testAnn.vcf

java -Xmx4g -jar /data/biocore/programs/snpEff/snpEff.jar \
  -v Pseudomonas_aeruginosa_pao1 /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1/pao1.genotyping.genotyped.gatk.g.vcf \
  > pao1.filtered.testAnn.vcf
```

Then, rename the VCF files where necessary.

```
cd /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1
for i in $(find ./ -type f -name "*.genotyped*.vcf"); do
  sed -e 's/NC_002516.2/Chromosome/g' "$i" > "${i%.vcf}".renamed.vcf
done
mkdir -p /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1/renamed
mv *renamed.vcf /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1/renamed

cd /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073
for i in $(find ./ -type f -name "*.genotyped*.vcf"); do
  sed -e 's/NZ_CP051263.1/Chromosome/g' "$i" > "${i%.vcf}".renamed.vcf
done
mkdir -p /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073/renamed
mv *renamed.vcf /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073/renamed

```

Finally, run snpEff itself for the annotations.

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/genotyping-pao1 \
  -l "pao1.genotyping.genotyped.filtered.gatk.g.renamed.vcf pao1.genotyping.genotyped.gatk.g.renamed.vcf" \
  -r Pseudomonas_aeruginosa_pao1 \
  -o /data/gencore/analysis_projects/8040442_Yang/snpEff-genotyping-pao1 \
  -s /data/biocore/programs/snpEff

sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/genotyping-cft073 \
  -l "cft073.genotyping.genotyped.gatk.g.renamed.vcf cft073.genotyping.genotyped.filtered.gatk.g.renamed.vcf" \
  -r Escherichia_coli_cft073 \
  -o /data/gencore/analysis_projects/8040442_Yang/snpEff-genotyping-cft073 \
  -s /data/biocore/programs/snpEff
```

### Individual VCF Annotation

I don't know how useful this will actually be, since I think the genotyping information will be more valuable for this project, but it could potentially be used as back-up or confirmatory information regarding specific timepoints.

First, identify the chromosome names used by the snpEff databases to ensure they match the chromosome names used in the output vcf files:

```
java -Xmx4g -jar /data/biocore/programs/snpEff/snpEff.jar \
  -v Pseudomonas_aeruginosa_pao1 /data/gencore/analysis_projects/8040442_Yang/variants-pao1/CAN1103.filtered.vcf \
  > CAN1103.filtered.testAnn.vcf
```

Then, rename the VCF files where necessary.

```
cd /data/gencore/analysis_projects/8040442_Yang/variants-pao1
for i in $(find ./ -type f -name "*.vcf"); do
  sed -e 's/NC_002516.2/Chromosome/g' "$i" > "${i%.vcf}".renamed.vcf
done
mkdir -p /data/gencore/analysis_projects/8040442_Yang/variants-pao1/renamed
mv *renamed.vcf /data/gencore/analysis_projects/8040442_Yang/variants-pao1/renamed

cd /data/gencore/analysis_projects/8040442_Yang/variants-cft073
for i in $(find ./ -type f -name "*.vcf"); do
  sed -e 's/NZ_CP051263.1/Chromosome/g' "$i" > "${i%.vcf}".renamed.vcf
done
mkdir -p /data/gencore/analysis_projects/8040442_Yang/variants-cft073/renamed
mv *renamed.vcf /data/gencore/analysis_projects/8040442_Yang/variants-cft073/renamed

```

Finally, run snpEff itself for the annotations.

```
sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/variants-pao1/renamed \
  -r Pseudomonas_aeruginosa_pao1 \
  -o /data/gencore/analysis_projects/8040442_Yang/snpEff-variants-pao1 \
  -s /data/biocore/programs/snpEff

sbatch /data/gencore/analysis_projects/8040442_Yang/scripts/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/8040442_Yang/variants-cft073/renamed \
  -r Escherichia_coli_cft073 \
  -o /data/gencore/analysis_projects/8040442_Yang/snpEff-variants-cft073 \
  -s /data/biocore/programs/snpEff
```
