## Standard Differential Expression Analysis

This analysis was done, unless otherwise specified, using the code and parameters found in the Bioinformatics Core's RNA-DEGs-Modular pipeline, version 1.0.1 released 03-14-2024, on ASU's SOL supercomputer.

### Alignment to Human Genome

Alignment of fastq reads to the human genome was carried out with STAR, first trying the 'default' parameters (found in the supplemental files folder for this GitHub release).

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A/A.alignment-wrapper.sh \
  --fastqDir /data/gencore/analysis_projects/7232055_Nolz/fastq \
  --intermediateDir /scratch/kawoodbu \
  --starParams /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/supplemental_files/star-params-default.txt \
  --refDir /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes \
  --alignmentOutDir /data/gencore/analysis_projects/7232055_Nolz/align-default \
  --quantOutDir /data/gencore/analysis_projects/7232055_Nolz/stringtie-default \
  --scriptDir /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A
  ```

With the default parameters, a substantial fraction of the reads were unmapped because the alignments were too short (< 66% of the total read length). So, alignment was rerun using the 'short' parameters, which set the minimum alignment percentage to only 30% of the total read length. The mismatch rate per base was near the range of normal RNA seq (0.57-0.71%, when the author of STAR notices a typical range of 0.3-0.6%), so retaining the typical max mismatch settings shouldn't reduce the quality of the output.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A/A.alignment-wrapper.sh \
  --fastqDir /data/gencore/analysis_projects/7232055_Nolz/fastq \
  --intermediateDir /scratch/kawoodbu \
  --starParams /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/supplemental_files/star-params-short.txt \
  --refDir /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes \
  --alignmentOutDir /data/gencore/analysis_projects/7232055_Nolz/align-short \
  --quantOutDir /data/gencore/analysis_projects/7232055_Nolz/stringtie-short \
  --scriptDir /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A
```

While alignment statistics did improve using the 'short' parameters, several of the samples still had very high '% of reads unmapped: too short' (up to 70%!). So I also ran the alignment using the 'lowqual' parameters, which require an alignment length of at least 20bp but doesn't take into consideration the ratio of alignment length to read length.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A/A.alignment-wrapper.sh \
  --fastqDir /data/gencore/analysis_projects/7232055_Nolz/fastq \
  --intermediateDir /scratch/kawoodbu \
  --starParams /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/supplemental_files/star-params-lowqual.txt \
  --refDir /data/gencore/databases/reference_genomes/human/Homo_sapiens_GRCh38_102/star_2.7.10a_indexes \
  --alignmentOutDir /data/gencore/analysis_projects/7232055_Nolz/align-lowqual \
  --quantOutDir /data/gencore/analysis_projects/7232055_Nolz/stringtie-lowqual \
  --scriptDir /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1/Module_A
```

Using the 'lowqual' parameters allowed for a much larger percentage of reads to be uniquely mapped, which combined with the multimapped reads (less than 10 mappings) meant that almost all reads were mapped in some usable way.

![Image](./markdown-images/STAR_QC_Reads.png)

### Differential Gene Expression

For basic differential expression based on the gene count matrix, I used the core's wrapper script for DEseq2, edgeR, and NOISeq analysis. This requires a comparisons file delineating which samples are part of each group and how the groups should be compared to each other.

```
bash /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.1-branch/Module_B/B1.DEG.rscripts.sh \
        -d /data/gencore/analysis_projects/7232055_Nolz/DEGS-branch \
        -g /data/gencore/analysis_projects/7232055_Nolz/stringtie-lowqual/gene_count_matrix.csv \
        -c /data/gencore/analysis_projects/7232055_Nolz/comparisons.csv \
        -s "deseq2 edger noiseq"
```

I noticed that this automatic pipeline was producing strange or no results, so I went through all the DEG R scripts individually. I adjusted to allow for single-column comparisons file (this was a major problem for this project) as well as a tryCatch clause to ensure that NOISeq was using the best normalization method for the available data. Since I was already running the scripts manually, I also assigned custom colors to the MDS plot in R to give that figure better readability.

To merge the genes identified with each tool and create Venn diagrams of the overlap, I used our script ``merge_de_both.py``:

```
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py ../comparisons.csv 0
```

NOISeq found no differential genes, edgeR found one for each condition, and DESeq2 found around 10-20 for each condition. I'm taking all into functional enrichment but it may not be enough to identify any functions of interest.
