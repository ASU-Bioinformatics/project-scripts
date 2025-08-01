# Code Used for Differential Gene Expression Analysis for iLab 8363170

## Module A: Alignment and Quantification

These samples are dual culture E. coli and Pseudomonas aeruginosa. I don't have STAR indexes made for either species, so I'll be running genome generate twice.
For consistency with past projects from Jiseon, I'm using the CFT073 E. coli strain instead of the K-12 standard.

The first attempt at building the E. coli indexes failed because there was a space in column two for most entries, which threw off STAR's ability to detect which features were CDS. I used the code `sed -i 's/Protein[[:space:]]Homology/ProteinHomology/g' CFT073.genomic.gtf` to remove that extra space.

```
# escherichia coli

refDir="/data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2"
fasta="$refDir"/"CFT073.genomic.fna"
gtf="$refDir"/"CFT073.genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon CDS \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 3000000000000
```

```
# pseudomonas aeruginosa

refDir="/data/gencore/databases/reference_genomes/pseudomonas_aeruginosa"
fasta="$refDir"/"GCF_000006765.1_PAO1.genomic.fna"
gtf="$refDir"/"GCF_000006765.1_PAO1.genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon CDS \
  --limitGenomeGenerateRAM 3000000000000
```
The naming for these samples is non-traditional so I'm renaming them.

```
cd /data/gencore/analysis_projects/8363170_Yang/fastq/
rename '_L005' '_SRN_L001' /data/gencore/analysis_projects/8363170_Yang/fastq/*
```

I'm using the GitHub repository for the RNA scripts this time, now that it's connected to Sol!

```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-default.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-default \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-default \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-default.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-default \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-default \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

Unfortunately, most reads showed up as too short or as multimapped for this culture, so I'm going to run a different parameter setting to try to capture those.


```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_ecoli \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

Again here most of the reads were multimapped. STAR will use the multimapped reads, but I'm concerned based on the percentages that many of these reads are binding multiple times to places on the genome that are very similar between species. For this reason I'm going to require slightly longer binding length but reduce the multimapping max number.

```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_ecoli \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-short30 \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30 \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa
# I forgot to change the output folder name for this one so it overwrote the regular bacterial short parameters
# but the results from that were so bad I don't think it really matters

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

There was an unexplainable error for sample Spx21-Ground-Mid-G2AG in the Pseudomonas alignment (it said biocore-rna wasn't a conda environment, when it had worked for all other samples). So, I have to rerun that sample separately and then restart the merge and quant scripts.

```
sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A1.cpu-align.sh \
  -s Spx21-Ground-Mid-G2AG \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -o /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short

for i in $(find /data/gencore/analysis_projects/8363170_Yang/fastq -type f -name "*fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq)
do
  echo $i
  sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A2.stringtie-quant.sh \
    -s "$i" \
    -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short30 \
    -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short30 \
    -g /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.genomic.gtf
done

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short30 \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short30
```

## Module B: Differential Expression

```
# e. coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8363170_Yang/ecoli-differentials-bacterial-short30 \
  -g /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8363170_Yang/comparisons.csv \
  -s "deseq2 edger noiseq"

# p. aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8363170_Yang/ecoli-differentials-bacterial-short30 \
  -g /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8363170_Yang/comparisons.csv \
  -s "deseq2 edger noiseq"

```
