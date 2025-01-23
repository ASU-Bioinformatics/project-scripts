# Module K
This module is a collection of (mostly) R scripts for differential expression analysis given the transcript count matrix from the abundance profiling in Module I. Because they can demand a lot of memory, especially edgeR, a wrapper script is included from which any of the included differential R scripts can be launched to the slurm scheduler on Sol. The resulting DE gene lists are merged, annotated, and filtered so that only DE genes identified by at least two tools are included. KEGG terms associated with the DE genes are extracted.

At this point, the only functional profile I've found that can handle KEGG terms across the multiple species represented in a microbiome is MicrobiomeProfiler, which functions largely as a Shiny app. Using the extracted KEGG terms, the Shiny app can be used to create barplots and dotplots of significant KEGG modules for each comparison, while a standard R function can be used to extract all tabular data corresponding to the KEGG significance analysis.

## K1.rscripts.wrapper.sh

Script K1 allows all the differential expression R scripts to be run through the Slurm scheduler on Sol, because they can take some time and require more memory than a standard RNA sequencing analysis. In addition to the input files needed for each script, the script name needs to be submitted to the wrapper.

Several of these packages are not currently publicly available on Sol, so they will need to be installed in each user's personal R library.

R packages needed for all scripts:
* `dplyr`
* `optparse`
* `viridis`


Input parameters:
* `-s,--script`: the name of the R script that should be called. Acceptable options include deseq2, edgeR, NOIseq, and clusterProfiler, and the script name is case-insensitive.
* `-d,--directory`: the name of the directory to which output files should be written. If an absolute pathname is not provided for the gene count matrix or comparisons file, the script will look for them in this folder.
* `-g,--geneMatrix`: the name of the file (absolute or within `--directory`) that contains the gene count matrix.
* `-c,--comparisons`: the name of the file (absolute or within `--directory`) that contains the pairwise comparison data for differential analysis. A sample comparison file is included here; essentially, sample names should be row names and comparison names should be the column headers. Within each comparison, each sample in the control group is represented by a 1 and each sample in the experimental group is represented by a -1. Any unused samples for a given comparison are represented as 0. The same file format can be used for all three differential expression scripts.

### K2.deseq2.R

This script uses the R package [DEseq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8,"Title") to normalize counts and determine differentially expressed genes using logistic regression models to calculate expression change and the Wald test to establish significance. If you choose to run only a single differential expression tool, this is the best one for use cases where each group has at least two samples.

Additional R packages needed:
* `DESeq2`
* `tibble`
* `ggplot2`
* `RColorBrewer`
* `ggrepel`
* `pheatmap`
* `stringr`

An example of how to call this script from the wrapper:

    sbatch ./I5.rscripts.wrapper.sh --script deseq2 \
      --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression \
      --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/gene.count.matrix.tsv \
      --comparisons /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/comparisons.csv

### K3.edgeR.R

This script uses the R package [edgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458?login=true,"Title"), designed for low-replicate datasets, to normalize counts and determine differentially expressed genes. EdgeR models gene expression with a negative binomial distribution and assesses differential expression using an adapted Fisher's exact test. This script has the highest memory demand, so if time and memory are an issue it can be skipped.

Additional R packages needed:
* `edgeR`
* `pheatmap`

An example of how to call this script from the wrapper:

    sbatch ./K1.rscripts.wrapper.sh --script edger \
      --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression \
      --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/gene.count.matrix.tsv \
      --comparisons /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/comparisons.csv

### K4.noiseq.R

This script uses the R package [NOISeq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/,"Title") to normalize counts and determine differentially expressed genes. NOIseq filters low-count features taking the experimental design into consideration and corrects for batch effect as part of normalization. For data sets with groups containing only a single sample, this is the best script to use for differential analysis as it was designed for that use case. However, in a recent comparison of DE R tools, this one was not among the top scorers, and I'll hopefully be replacing it soon.

Additional R packages needed:
* `NOISeq`
* `tidyverse`

An example of how to call this script from the wrapper:

    sbatch ./K1.rscripts.wrapper.sh --script noiseq \
      --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression \
      --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/gene.count.matrix.tsv \
      --comparisons /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression/comparisons.csv

## K5.merge-annotate-degs.sh and K5.merge_de.py

This script takes the differential expression gene lists from scripts K2-K4 and merges them into summary gene lists containing data from however many tools were used. It will output a summary gene list, Venn diagram, and a stringent gene list (genes flagged by at least two tools) for each pairwise comparison.

In addition, it will annotate the gene lists and extract KEGG terms associated with the stringent gene lists for downstream functional profiling.

## K6.functionalProfiler.R

Unlike the other scripts in this module, this one is manual and fairly time-intensive, as the R package implemented here functions as a Shiny app from which figures need to be downloaded. The KEGG term lists for each comparison are tested for significance and dotplots and barplots are generated when significant KEGG modules are present. I typically require >2 KEGG terms to be present in the comparison to run the Shiny profiler.

Tabular output for the significance testing can be run with a loop in the standard R console, so is much faster and can be done for all comparisons (although comparisons with no KEGG terms identified from the significant gene lists will result in empty tables).

## Before moving on:

Look at the MDS/PCA plots and heatmaps to see if there is any clustering between the samples. If there isn't any clear clustering, it's unlikely that the sample groups are significantly different in their expressed functions, despite the presence of individual differentially expressed genes.

If the two primary components of the MDS/PCA plots don't convey a majority of the difference between the samples, it's even less likely that there are meaningful differences between sample groups. This can mean that all we're seeing is the normal biological variation between individual samples. When this is the case, the heatmaps should reflect this, with high or low expression patterns differing by sample rather than showing any group similarities.

Functional analysis can show the most interesting results when functional differences are present; however, any comparisons including a group with only one sample should be treated with some doubt.
