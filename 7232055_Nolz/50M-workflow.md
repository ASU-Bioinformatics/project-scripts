## Standard Differential Expression Analysis

This analysis was done, unless otherwise specified, using the code and parameters found in the Bioinformatics Core's RNA-DEGs-Modular pipeline, version 1.0.1 released 03-14-2024, on ASU's SOL supercomputer.

### Alignment to Human Genome

Alignment of fastq reads to the human genome was carried out with STAR. Based on the mapping results from the preliminary data, I also ran the alignment using the 'lowqual' parameters, which require an alignment length of at least 20bp but doesn't take into consideration the ratio of alignment length to read length.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  --fastqDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/fastq \
  --intermediateDir /scratch/kawoodbu \
  --starParams /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-lowqual.txt \
  --refDir /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes \
  --alignmentOutDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/align-lowqual \
  --quantOutDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual \
  --scriptDir /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A
```

Due to an error with the fastq.gz file for 10-39, I needed to rerun the alignment and quant for that sample only, and then rerun the merge quant.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A1.cpu-align.sh \
  --sampleID 10-39 --intermediateDir /scratch/kawoodbu \
  --fastqDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/fastq \
  --starParams /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-lowqual.txt \
  --refDir /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes \
  --outputDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/align-lowqual

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A2.stringtie-quant.sh \
  --sampleID 10-39 \
  --alignmentDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/align-lowqual \
  --quantDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual \
  --refGTF /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes/Homo_sapiens.GRCh38.102.gtf

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A3.merge-quant.sh \
  --scriptDir /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A \
  --quantDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual \
  --alignmentDir /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/align-lowqual
```

![Image](./markdown-images/STAR_QC_Reads.png)

### Differential Gene Expression

For basic differential expression based on the gene count matrix for the MGX samples, I had used the core's wrapper script for DEseq2, edgeR, and NOISeq analysis. This requires a comparisons file delineating which samples are part of each group and how the groups should be compared to each other.

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
        -d /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/DEGS-lowqual \
        -g /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual/gene_count_matrix.csv \
        -c /data/gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/comparisons.csv \
        -s "deseq2 edger noiseq"
```

I noticed that this automatic pipeline was producing strange or no results, so I went through all the DEG R scripts individually. I adjusted to allow for single-column comparisons file (this was a major problem for this project) as well as a tryCatch clause to ensure that NOISeq was using the best normalization method for the available data. Since I was already running the scripts manually, I also assigned custom colors to the MDS plot in R to give that figure better readability.

Instead of trying to run the automatic pipeline for the MTX samples, I just ran them manually to begin with.

To merge the genes identified with each tool and create Venn diagrams of the overlap, I used our script ``merge_de_both.py``:

```
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py ../comparisons.csv 0
```

NOISeq found no differential genes, edgeR found one for each condition, and DESeq2 found around 10-20 for each condition. I'm taking all into functional enrichment but it may not be enough to identify any functions of interest.
