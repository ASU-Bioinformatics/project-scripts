# RNA Sequencing
## Explanation of Methods and Output
### Folder Structure

```
-Project ID
|--project_id_overview.pdf
|--1.fastq
  |--sample1_S01_L001_R1_001.fastq.gz
  |--sample1_S01_L001_R2_001.fastq.gz
  |--etc.
|--2.sequencing-qc
  |--multiqc_report.html
  |--fastqc
    |--sample1_S01_L001_R1_001_fastqc.html
    |--sample1_S01_L001_R1_001_fastqc.zip
    |--sample1_S01_L001_R2_001_fastqc.html
    |--sample1_S01_L001_R2_001_fastqc.zip
    |--etc.
  |--multiqc_data
    |--multiqc_citations.txt
    |--multiqc_data.json
    |--multiqc_fastqc.txt
    |--multiqc_general_stats.txt
    |--multiqc_software_versions.txt
    |--multiqc_sources.txt
    |--multiqc.log
|--3.alignment-and-counts
  |--CPM.txt
  |--deseq2_normalized_counts.txt
  |--gene_count_matrix.csv
  |--gene_fpkm.tsv
  |--gene_tpm.tsv
  |--STAR_QC_Errors.png
  |--STAR_QC_Reads.png
  |--STAR_Stats.csv
  |--transcript_count_matrix.csv
  |--transcript_fpkm.tsv
  |--transcript_tpm.tsv
  |--bams
    |--sample1_STARAligned.sortedByCoord.out.bam
    |--etc.
|--4.differential-expression
  |--allDEG.heatmap.pdf
  |--MDSplot.edgeR.pdf
  |--PCAplot.deseq.pdf
  |--PCAplot.noiseq.pdf
  |--top50genes.edgeR.heatmap.pdf
  |--top50genes.deseq.heatmap.pdf
  |--clustering
    |--clustered.heatmap.pdf
    |--clusters-by-samples.pdf
    |--genes-in-cluster.txt
    |--sumSquares.pdf
  |--group1_vs_group2 (one folder for each group comparison)
    |--group1_vs_group2.deg.down.venn.minLFC0.png
    |--group1_vs_group2.deg.up.venn.minLFC0.png
    |--group1_vs_group2.deg.down.venn.minLFC1.png
    |--group1_vs_group2.deg.up.venn.minLFC1.png
    |--group1_vs_group2.deg.down.venn.minLFC2.png
    |--group1_vs_group2.deg.up.venn.minLFC2.png
    |--group1_vs_group2.deseq.MAplot.pdf
    |--group1_vs_group2.deseq.volcano.png
    |--group1_vs_group2.edgeR.MAplot.pdf
    |--group1_vs_group2.noiseq.explot.pdf
    |--group1_vs_group2.noiseq.MDplot.pdf
    |--deseq-lists
      |--group1_vs_group2.deg.all.deseq.txt
      |--group1_vs_group2.deg.down.deseq.txt
      |--group1_vs_group2.deg.up.deseq.txt
    |--edgeR-lists
      |--group1_vs_group2.deg.all.edgeR.txt
      |--group1_vs_group2.deg.down.edgeR.txt
      |--group1_vs_group2.deg.up.edgeR.txt
    |--noiseq-lists
      |--group1_vs_group2.deg.all.noiseq.txt
      |--group1_vs_group2.deg.down.noiseq.txt
      |--group1_vs_group2.deg.up.noiseq.txt
    |--merged-deglists
      |--group1_vs_group2.deg.down.minLFC0.csv
      |--group1_vs_group2.deg.up.minLFC0.csv
      |--group1_vs_group2.deg.down.minLFC1.csv
      |--group1_vs_group2.deg.up.minLFC1.csv
      |--group1_vs_group2.deg.down.minLFC2.csv
      |--group1_vs_group2.deg.up.minLFC2.csv
|--5.functional-enrichment
  |--group1_vs_group2 (one folder for each group comparison)
    |--go-enrichment
      |--go.cp.down.barplot.group1_vs_group2.pdf
      |--go.cp.down.cnetplot.group1_vs_group2.pdf
      |--go.cp.down.dotplot.group1_vs_group2.pdf
      |--go.cp.down.emapplot.group1_vs_group2.pdf
      |--go.cp.up.barplot.group1_vs_group2.pdf
      |--go.cp.up.cnetplot.group1_vs_group2.pdf
      |--go.cp.up.dotplot.group1_vs_group2.pdf
      |--go.cp.up.emapplot.group1_vs_group2.pdf
      |--go.group1_vs_group2.down.sig.csv
      |--go.group1_vs_group2.up.sig.csv
    |--kegg-pathways (only for organisms with KEGG pathway data available)
      |--hsa#####.group1_vs_group2.png (one for each significantly enriched pathway)
      |--kegg.group1_vs_group2.clusterProfiler.sig.csv
```
### Sequencing Reads and Quality Metrics

*Data found in folders `1.fastq` and `2.sequencing.qc`*

The subfolder `1.fastq` contains the demultiplexed fastq reads from the Illumina sequencer, including the forward and reverse reads as separate files. `2.sequencing.qc` contains the fastqc quality metrics for each fastq file as well as a quality summary generated with [MultiQC]("https://multiqc.info/", "Title") that includes an HTML visualization and tabular statistical information.

### Alignment of Short Read Sequences

*Data found in folder `3.alignment-and-counts`*

The first step of RNA-sequencing analysis is to align the short reads to a reference genome to capture read counts for transcripts/genes despite their variable lengths. If a reference genome is available in NCBI, we can use this as the alignment reference; otherwise, we can use a personal genome assembly from the originating researcher. In either case, to include annotations for each gene and enable functional enrichment analysis downstream, both a fasta and gff/gtf file are required.

The alignment itself is performed with the Alex Dobin's open-source tool [STAR]("https://github.com/alexdobin/STAR", "Title"), which searches for the maximum mappable portion of each read on the genome iteratively for each subsequent unmapped portion, enabling it to detect splice junctions. After the maximum mappable seeds are found, they are stitched together to incorporate alignments with mismatches, indels, and splice junctions into full-length transcripts.

The STAR analysis outputs a BAM alignment file for each sample. BAM files are binary computer-readable files that contain all alignment information for query-subject pair identified by the aligner (in this case, STAR). They are primarily used as intermediate files to feed into the next step of the pathway, but can be converted to human-readable SAM files if desired for manual inspection.

Once the alignments are generated, we use the open-source tool [StringTie]("https://ccb.jhu.edu/software/stringtie/", "Title"), with several publicly-available and in-house helper scripts, to count all transcripts and genes based on the reference genome annotations and generate a unified gene count matrix (the accompanying transcript count matrix can be used for determining alternate splice junctions and variants for individual genes, but that analysis is not included in out standard RNA-sequencing deliverables).

The output from StringTie is provided in the folder `3.alignment-and-counts` and includes the transcript and gene count matrices as well as FPKM and TMM normalized copies of those matrices.

### Differential Expression Analysis

*Data primarily found in folder `4.differential-expression`*

Once the gene count matrix has been generated, it is normalized and analyzed for differential expression. Three primary open-source R packages ([DEseq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8,"Title"), [edgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458?login=true,"Title"), and [NOISeq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/,"Title")) are used for this process by the core, as all have different strengths and weaknesses and approach the statistical problems in slightly different ways. DEseq2 determines differentially expressed genes using logistic regression models to calculate expression change and the Wald test to establish significance. EdgeR is designed for low-replicate datasets (less than 7-8 replicates per condition); it models gene expression with a negative binomial distribution and assesses differential expression using an adapted Fisher's exact test. NOIseq filters low-count features taking the experimental design into consideration and corrects for batch effect as part of normalization; it was originally designed for experiments with no replicates and is ideal for that use case.

DEseq2's normalization algorithm is considered one of the top metrics for RNA sequencing, because it takes into consideration both read depth (like FPKM, CPM, and TPM metrics) and population composition (the change in gene expression ratios due to increased reads for some genes). The resulting normalized gene count matrix is included in the folder `3.alignment-and-counts` along with the CPM-normalized matrix created by edgeR (downstream edgeR analysis uses a similar normalization method to DEseq2, not the CPM matrix, however). The remaining output files from these three programs are provided in the folder `4.differential-expression`.

For each tool, a PCA/MDS plot is included showing clustering patterns across the top two axes of variability. A heatmap representing the top 50 most variable genes (according to each tool's normalization algorithm) is also included, along with a heatmap of all differentially expressed genes.

For the highest certainty of results, select only those genes considered differential by all three tools. The core's standard process is to take only genes considered differentially expressed by at least two tools into clustering and functional enrichment analysis.

Within the `clustering` subfolder of `4.differential-expression`, there are the results of gene expression clustering analysis. DEseq2 normalized count values are extracted for genes determined to be differential by at least two of the R packages and grouped by similar patterns of expression across all samples. A heatmap showing cluster assignment, a scatter plot of average expression of each cluster for each sample, and a table containing the list of genes (by id) within each cluster are included.

Each group comparison has its own subfolder within `4.differential-expression`. These folders contain the differential gene lists generated by all three DEG analysis methods, the merged differential gene lists containing information for each gene from each tool (filtered at a minimum log2 fold change of 0, 1, or 2 respectively), Venn diagrams showing the overlap between differentially expressed genes determined by each tool, an MA plot and volcano plot from DEseq2, an MA plot from edgeR, and an MD plot and expression plot from NOIseq.

(A brief explanation of the plot types:
  * MA plot: these plots graph the average log2 fold change between sample groups for each gene against the average normalized count for each gene, showing the relationship between variability and expression level for each gene in a given group comparison.
  * Volcano plot: these plots show the general relationship between log2 fold change and adjusted p-value (graphed as -log10 of the value so that increased significance is higher on the y-axis). They typically look like an erupting volcano because the significance of a potential DEG call tends to increase as the calculated log2 fold change in expression increases.
  * MD plot: this plot graphs the additive difference between gene expression for each gene in the two groups against the multiplicative difference between for that same comparison. This provides a visual of how differential expression tracks across expression level of each gene.
  * Expression plot: this plot graphs the average expression of each gene in group 2 against the average expression of each gene in group 1, and typically looks fairly linear with significant outliers marked in red.

  )

### Functional Enrichment Analysis

*Data found in folder `5.functional-enrichment`*

Each pairwise comparison between two sample groups has its own subfolder of results. If an organism has no available functional annotations (for example, a de novo assembly), this entire folder will be missing. If an organism has only GO annotations without KEGG annotations (most less-common organisms), only the files from GO enrichment will be included. For the 1,309 eukaryotes and 8,898 prokaryotic genomes in the KEGG database, a KEGG pathway enrichment analysis will also be carried out.

Both GO and KEGG analyses are performed with the R package [clusterProfiler]("https://guangchuangyu.github.io/software/clusterProfiler/" "Title"). For GO enrichment, returned files include the lists of significant GO terms identified in the data (simplified into representative terms, since higher-level GO terms can be statistically significant just because one or two of their children are over-represented) as well as bar plots, network plots that link differentially expressed genes with differentially expressed GO terms, mapping plots that join GO term nodes based on the similarity of gene expression patterns, and dot plots for over and under-expressed genes for each comparison.

For KEGG enrichment, an image of each significantly enriched pathway is provided (with differentially expressed genes highlighted) as well as a CSV file containing metrics for all significantly enriched pathways.

## Additional Information

For more information or to schedule a consultation with the core, please email us at <asubioinformatics@asu.edu>. If you'd like additional or customized analyses performed on your data, please contact us at the same email for an updated quote.
