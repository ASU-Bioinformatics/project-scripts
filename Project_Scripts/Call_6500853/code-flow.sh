sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/A1.alignment.sh \
  --refDir /data/gencore/databases/reference_genomes/drosophila_melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT \
  --refID GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta \
  --fastqDir /data/gencore/analysis_projects/6500853_Call/fastq \
  --alignmentDir /data/gencore/analysis_projects/6500853_Call/alignment

java -jar picard.jar CreateSequenceDictionary \
      R=/data/gencore/databases/reference_genomes/drosophila_melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta \
      O=/data/gencore/databases/reference_genomes/drosophila_melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.dict

sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/A2.variant-calling.sh \
  --inputDir /data/gencore/analysis_projects/6500853_Call/alignment \
  --ref /data/gencore/databases/reference_genomes/drosophila_melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta \
  --outDir /data/gencore/analysis_projects/6500853_Call/variants

# edit chromosome names and remove extraneous scaffolds for Drosophila alignment that I messed up
vcftools --vcf hf_de_addrg_Park-Homozygous-Male2_bwa_aligned.sorted.vcf \
         --out hf_de_addrg_Park-Homozygous-Male2_bwa_aligned.sorted.cleaned \
         --chr NT_033779.5 --chr NT_033778.4 --chr NT_037436.4 --chr NT_033777.3 \
         --chr NC_004353.4 --chr NC_024512.1 --chr NC_004354.4 --recode

grep -v '^##contig=<ID=NW_' hf_de_addrg_Park-Homozygous-Male2_bwa_aligned.sorted.cleaned.recode.vcf > hf_de_addrg_Park-Homozygous-Male2_bwa_aligned.sorted.cleaned.vcf


# after this I renamed the alignment and variant files to match the new naming conventions

i=Park-Homozygous-Male2
sed -i 's/NC_024511.2/chrM/g' "$i".filtered.vcf
sed -i 's/NC_004354.4/chrX/g' "$i".filtered.vcf
sed -i 's/NC_004353.4/chr4/g' "$i".filtered.vcf
sed -i 's/NC_024512.1/chrY/g' "$i".filtered.vcf
sed -i 's/NT_033779.5/chr2L/g' "$i".filtered.vcf
sed -i 's/NT_033778.4/chr2R/g' "$i".filtered.vcf
sed -i 's/NT_037436.4/chr3L/g' "$i".filtered.vcf
sed -i 's/NT_033777.3/chr3R/g' "$i".filtered.vcf

grep -v '^##contig=<ID=NW_' "$i".filtered.vcf > "$i".filtered1.vcf

vcftools --vcf "$i".filtered1.vcf \
         --out "$i".filtered2.vcf \
         --chr chrX --chr chr4 --chr chrY --chr chr2L \
         --chr chr2R --chr chr3L --chr chr3R --recode

bedtools genomecov -ibam Park-Homozygous-Male2.sorted.rdgrp.mrkd.bam -max 10 > Park-Homozygous-Male2.genomecov.txt
grep -v '^NW_' Park-Homozygous-Male2.genomecov.txt > Park-Homozygous-Male2.cleaned.genomecov.txt
grep -v 'NC_024511.2' Park-Homozygous-Male2.cleaned.genomecov.txt > Park-Homozygous-Male2.renamed.genomecov.txt

sed -i 's/NC_004354.4/chrX/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NC_004353.4/chr4/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NC_024512.1/chrY/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NT_033779.5/chr2L/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NT_033778.4/chr2R/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NT_037436.4/chr3L/g' Park-Homozygous-Male2.renamed.genomecov.txt
sed -i 's/NT_033777.3/chr3R/g' Park-Homozygous-Male2.renamed.genomecov.txt

i=Park-Homozygous-Male2
java -Xmx16g -jar /data/biocore/programs/snpEff/snpEff.jar \
                  Drosophila_melanogaster "$i".filtered2.vcf \
                  -s "$i".filtered.snpEff.html > "$i".filtered.snpEff.vcf

bcftools view -i 'ANN ~ "HIGH"' Park-Homozygous-Male2.filtered.snpEff.vcf > Park-Homozygous-Male2.filtered.snpEff.high-impact.vcf

grep 'ANN=' Park-Homozygous-Male2.filtered.snpEff.high-impact.vcf | awk '{print $6"|"$8}' | awk -F '|' '{print $1"\t" $6}' > high-impact-genes.txt

# running this code will output a four column tab-delimited file where:
# col1 = ensembl gene id
# col2 = reads for reference genome allele
# col3 = reads for mutated allele
# col4 = total reads for position
grep 'ANN=' Park-Homozygous-Male2.filtered.snpEff.high-impact.vcf | awk '{print $8"|"$10}' | awk -F '|' '{print $5":"$NF}' | awk -F ":" '{print $1"\t"$3"\t"$4}' | awk -F "," '{print $1"\t"$2}' > high-impact-gene-read-depths.txt

###
# ran gridss for large structural variant detection

vcftools --vcf Park-Homozygous-Male2.gridss.vcf \
         --out Park-Homozygous-Male2.gridss.cleaned.vcf \
         --chr NT_033779.5 --chr NT_033778.4 --chr NT_037436.4 \
         --chr NT_033777.3 --chr NC_004353.4 --chr NC_024512.1 \
         --chr NC_004354.4 --recode

grep -v '^##contig=<ID=NW_' Park-Homozygous-Male2.gridss.cleaned.vcf.recode.vcf > Park-Homozygous-Male2.gridss.cleaned.vcf

# after this I renamed the alignment and variant files to match the new naming conventions

sed -i 's/NC_024511.2/chrM/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NC_004354.4/chrX/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NC_004353.4/chr4/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NC_024512.1/chrY/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NT_033779.5/chr2L/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NT_033778.4/chr2R/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NT_037436.4/chr3L/g' Park-Homozygous-Male2.gridss.cleaned.vcf
sed -i 's/NT_033777.3/chr3R/g' Park-Homozygous-Male2.gridss.cleaned.vcf

# remove all calls that don't pass the built in quality filter

vcftools --vcf Park-Homozygous-Male2.gridss.cleaned.vcf \
         --out Park-Homozygous-Male2.chr3.gridss.hardfiltered.vcf \
         --remove-filtered-all --chr chr3L --chr chr3R --recode

# now rerun snpeff on the vcf file containing all the breakend information.
java -Xmx16g -jar /data/biocore/programs/snpEff/snpEff.jar \
                  Drosophila_melanogaster Park-Homozygous-Male2.chr3.gridss.hardfiltered.vcf \
                  -s Park-Homozygous-Male2.chr3.gridss.hardfiltered.snpEff.html > Park-Homozygous-Male2.chr3.gridss.hardfiltered.snpEff.vcf

# filter only high impact variants
bcftools view -i 'ANN ~ "HIGH"' Park-Homozygous-Male2.chr3.gridss.hardfiltered.snpEff.vcf > Park-Homozygous-Male2.chr3.gridss.hardfiltered.snpEff.high-impact.vcf

# ran delly for large structural variant detection
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/delly-env

delly call -g /data/gencore/databases/reference_genomes/drosophila_melanogaster/GCF_000001215.4_Release_6_plus_ISO1_MT/GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fasta \
           -o Park-Homozygous-Male2.bcf ./Park-Homozygous-Male2.sorted.rdgrp.mrkd.bam

module load bcftools-1.10.2-gcc-12.1.0
bcftools convert Park-Homozygous-Male2.bcf -o Park-Homozygous-Male2.vcf

module load vcftools-0.1.14-gcc-11.2.0
vcftools --vcf Park-Homozygous-Male2.vcf \
         --out Park-Homozygous-Male2.cleaned.vcf \
         --chr NT_033779.5 --chr NT_033778.4 --chr NT_037436.4 \
         --chr NT_033777.3 --chr NC_004353.4 --chr NC_024512.1 \
         --chr NC_004354.4 --recode

grep -v '^##contig=<ID=NW_' Park-Homozygous-Male2.cleaned.vcf.recode.vcf > Park-Homozygous-Male2.cleaned.vcf

# after this I renamed the alignment and variant files to match the new naming conventions

sed -i 's/NC_024511.2/chrM/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NC_004354.4/chrX/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NC_004353.4/chr4/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NC_024512.1/chrY/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NT_033779.5/chr2L/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NT_033778.4/chr2R/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NT_037436.4/chr3L/g' Park-Homozygous-Male2.cleaned.vcf
sed -i 's/NT_033777.3/chr3R/g' Park-Homozygous-Male2.cleaned.vcf

# remove all calls that don't pass the built in quality filter

vcftools --vcf Park-Homozygous-Male2.cleaned.vcf \
         --out Park-Homozygous-Male2.chr3.hardfiltered.vcf \
         --remove-filtered-all --chr chr3L --chr chr3R --recode

module load jdk-12.0.2_10-gcc-12.1.0
java -Xmx16g -jar /data/biocore/programs/snpEff/snpEff.jar \
                  Drosophila_melanogaster Park-Homozygous-Male2.chr3.hardfiltered.vcf.recode.vcf \
                  -s Park-Homozygous-Male2.chr3.hardfiltered.delly.snpEff.html > Park-Homozygous-Male2.chr3.hardfiltered.delly.snpEff.vcf
