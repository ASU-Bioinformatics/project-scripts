# Code Used for Project 7763550 Analysis

## Alignment to Hare Genome

### Create STAR Indexes for Hare/Myxoma Dual Reference

I used the current reference genome, GCF_033115175.1_mLepTim1.pri, as the source for the Hare STAR indexes. For myxoma, I used the same genome as for the previous human/myxoma dual reference (literally, I grabbed the intermediate files so I wouldn't have to remake them!)

```
cd /data/gencore/databases/reference_genomes/hare/hare-myxoma-dual
cat GCF_033115175.1_mLepTim1.pri_genomic.gtf myxoma-AF170726.2.gtf > hare-myxoma-dual.gtf
cat GCF_033115175.1_mLepTim1.pri_genomic.fna myxoma-AF170726.2.fasta > hare-myxoma-dual.fasta

#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.hare.out
#SBATCH -o slurm.%j.starGG.hare.err

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/hare/hare-myxoma-dual"
fasta="$refDir"/"hare-myxoma-dual.fasta"
gtf="$refDir"/"hare-myxoma-dual.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon exon \
  --limitGenomeGenerateRAM 3000000000000
```

My previous notes on creating a dual reference with myxoma imply that I'll have to manually add some annotations (I think to the geneInfo.tab file), but I should be able to use the intermediate files I created at that time.

### Create STAR Indexes for Hare Reference

I used the current reference genome, GCF_033115175.1_mLepTim1.pri, as the source for the Hare STAR indexes.

```
#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.hare.out
#SBATCH -o slurm.%j.starGG.hare.err

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/hare"
fasta="$refDir"/"GCF_033115175.1_mLepTim1.pri_genomic.fna"
gff="$refDir"/"GCF_033115175.1_mLepTim1.pri_genomic.gff"
gtf="$refDir"/"GCF_033115175.1_mLepTim1.pri_genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon exon \
  --limitGenomeGenerateRAM 3000000000000
```

### Format Annotation Files for Functional Analysis

Because hare is a relatively uncommon reference genome to use, it doesn't have a GO annotation / gene linkup reference easily accessible in R (ie, there is no Org.xx.db package for it). So, I'm having to manually create a list of gene-to-GO terms that I can use for clusterProfiler functional analysis.

```
grep 'GO' GCF_033115175.1-RS_2024_01_gene_ontology.gaf | awk '{print $5,$3}' | tr ' ' '\t' > hare-goterms.txt
```

### Dual Genome Workflow

#### Align to Dual Reference

I always try alignment with the default parameters first, but for Junior those always seem to align at <70% at best.

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7763550_Junior_RNA/fastq \
  -i /scratch/kawoodbu/STAR-align-dual \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-default.txt \
  -r /data/gencore/databases/reference_genomes/hare/hare-myxoma-dual-STAR2.7.10a \
  -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-dual-default \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-dual-default \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A
```

So, I additionally aligned with the short parameters (for the samples I checked, >30% of reads were rejected for having too-short alignments).

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7763550_Junior_RNA/fastq \
  -i /scratch/kawoodbu/STAR-align-dual-short \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-short.txt \
  -r /data/gencore/databases/reference_genomes/hare/hare-myxoma-dual-STAR2.7.10a \
  -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-dual-short \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-dual-short \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A
```

#### Stringtie Quantification

Typically this should happen automatically using the alignment wrapper script; however, there was an error with the GTF file. I took some advice from the internet and ran the following script on the GTF file, then replaced the referenced GTF file in the genomeParameters file with the updated GTF.

```
awk '$3 != "gene" ' hare-myxoma-dual.gtf > hare-myxoma-dual_no_genes.gtf
```

After that, I ran the stringtie quant script independently. Because the short parameters had much better alignment percentages, I only continued with stringtie quantification for that set of bam files.

```
for i in Laus1 Laus3 Laus4 Mock1 Mock2 Mock5 Toc159KO1 Toc159KO2 Toc159KO5 Toc1 Toc2 Toc3
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A2.stringtie-quant.sh \
    -s "$i" -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-dual-short \
    -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-dual-short \
    -g /data/gencore/databases/reference_genomes/hare/hare-myxoma-dual-STAR2.7.10a/hare-myxoma-dual_no_genes.gtf
done;
```

Once the individual quants were complete, I ran the merge script independently.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-dual-short
```

#### DEG Scripts in R

The differential expression comparison scripts can be run automatically with the wrapper script, although the output should be carefully checked to ensure success.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/7763550_Junior_RNA/degs-dual-short/ \
  -g /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-dual-short/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/7763550_Junior_RNA/comparisons.csv \
  -s "deseq2 edger noiseq"
```


### Hare Only Workflow

#### Alignment (Full Wrapper Attempt)

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7763550_Junior_RNA/fastq \
  -i /scratch/kawoodbu/STAR-align-hare \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-default.txt \
  -r /data/gencore/databases/reference_genomes/hare/hare-STAR2.7.10a \
  -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-hare-default \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-default \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A
```

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7763550_Junior_RNA/fastq \
  -i /scratch/kawoodbu/STAR-align-hare-short \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-short.txt \
  -r /data/gencore/databases/reference_genomes/hare/hare-STAR2.7.10a \
  -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-hare-short \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-short \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A
```

#### Standalone Stringtie Quantification

Typically this should happen automatically using the alignment wrapper script; however, there was an error with the GTF file. I took some advice from the internet and ran the following script on the GTF file, then replaced the referenced GTF file in the genomeParameters file with the updated GTF.

```
awk '$3 != "gene" ' GCF_033115175.1_mLepTim1.pri_genomic.gtf > GCF_033115175.1_mLepTim1.pri_genomic_nogenes.gtf
```

After that, I ran the stringtie quant script independently. Because the short parameters had much better alignment percentages, I only continued with stringtie quantification for that set of bam files.

```
for i in Laus1 Laus3 Laus4 Mock1 Mock2 Mock5 Toc159KO1 Toc159KO2 Toc159KO5 Toc1 Toc2 Toc3
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A2.stringtie-quant.sh \
    -s "$i" -a /data/gencore/analysis_projects/7763550_Junior_RNA/alignment-hare-short \
    -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-short \
    -g /data/gencore/databases/reference_genomes/hare/hare-STAR2.7.10a/GCF_033115175.1_mLepTim1.pri_genomic_nogenes.gtf
done;
```

Once the individual quants were complete, I ran the merge script independently.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A \
  -q /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-short
```

#### DEG Scripts in R

The differential expression comparison scripts can be run automatically with the wrapper script, although the output should be carefully checked to ensure success.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/7763550_Junior_RNA/degs-hare-short/ \
  -g /data/gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-short/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/7763550_Junior_RNA/comparisons.csv \
  -s "deseq2 edger noiseq"
```
