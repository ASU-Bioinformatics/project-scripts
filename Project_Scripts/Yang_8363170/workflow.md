# Code Used for Differential Gene Expression Analysis for iLab 8363170

## Module A: Alignment and Quantification

These samples are dual culture E. coli and Pseudomonas aeruginosa. I don't have STAR indexes made for either species, so I'll be running genome generate twice.
For consistency with past projects from Jiseon, I'm using the CFT073 E. coli strain instead of the K-12 standard.

```
# escherichia coli

refDir="/data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1"
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
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1 \
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
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```
