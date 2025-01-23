# RNA Sequencing
This workflow can be run via a single command for the whole pipeline, separate commands for the alignment and analysis modules, or individual commands for each step of the process, providing flexibility for rerunning samples or testing parameters.

## Alignment Workflow

This workflow combines all individual scripts in module A, and includes alignment to the reference with [STAR]("https://github.com/alexdobin/STAR", "Title"), quantification of read counts per gene, and QC statistics for the alignment.

The set up for the workflow assumes the Sol environment at ASU, with a future release I'll include a better way to update sbatch parameters for other environments. I do have a yml file to install the mamba environment used by these scripts, although the name of the environment would need to be manually updated in your copy of the scripts.

#### Required Tools
* [STAR]("https://github.com/alexdobin/STAR", "Title") - on Sol, must be locally installed
* [Samtools]("https://www.htslib.org/", "Title") - on Sol, can be uploaded with ``module load samtools-1.16-gcc-11.2.0``
* [StringTie]("https://ccb.jhu.edu/software/stringtie/", "Title") and other python packages contained in biocoreRNA.yml - on Sol, can be accessed in the environment ``/data/biocore/programs/mamba-envs/biocore-rna``

#### Usage
    sbatch A.alignment-wrapper.sh \
      --fastqDir /path/to/fastqs \
      --intermediateDir /path/to/FIFO/directory \
      --starParams /path/to/star.parameters.txt \
      --refDir /path/to/reference/directory \
      --alignmentOutDir /path/to/alignment/output/directory \
      --quantOutDir /path/to/stringtie/output/directory \
      --scriptDir /path/to/auxiliary-script/directory \
      --list "file1 file2 file3 ..." \
      (--readTypeSingle)

#### Arguments
* ``-f,--fastqDir``: directory containing gzipped fastq files, named with standard Illumina formatting (ie, sid-1_S01_L001_R1.fastq.gz)
* ``-i,--intermediateDir``: FIFO-enabled server for STAR alignment
* ``-p,--starParams``: text file containing STAR parameters (see examples in scripts folder)
* ``-r,--refDir``: directory containing the reference genome, GTF file, and all STAR index files
* ``-a,--alignmentOutDir``: location for output alignment files (eg, sorted bam files)
* ``-q,--quantOutDir``: location for output stringtie quantification files (eg, transcript abundance tables for each sample)
* ``-s,--scriptDir``: directory where the auxiliary python and perl scripts can be found (scripts A4-A7)
* ``-t,--readTypeSingle``: optional; specifies that input reads are single-end (unidirectional)
* ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow
* ``-h,--help``: prints an informational message and exits script

## Individual Modules

### STAR Alignment on CPU Server

The script for this module is A1.cpu-align.sh, and relies on the [STAR]("https://github.com/alexdobin/STAR", "Title") alignment software by Alex Dobin and the [Samtools]("https://www.htslib.org/", "Title") suite for manipulating high-throughput sequencing data. For biocore use, STAR needs to be locally installed and Samtools can be accessed through the Mamba environment at ``/data/biocore/programs/mamba-envs/biocore-rna``, which is hardcoded into the script. If you don't have access to that environment, Samtools can also be loaded as a Sol module with the command ``module load samtools-1.16-gcc-11.2.0``.

When running this script independently of the alignment workflow, only one sample can be provided at a time. Alignment is time-intensive, so it is far more efficient to run them separately. You can of course submit the script within a for loop to accelerate the job submission if you want to use the script outside of the workflow wrapper!

#### Usage:
    sbatch A1.cpu-align.sh \  
      --sampleID "Massilia-plus-A" \
      --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
      --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
      --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-default.txt" \
      --refDir "/data/gencore/databases/reference_genomes/massilia" \
      --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment"``

#### Arguments:
* ``-s,--sampleID``: sample ID includes all text before the first underscore in a standard Illumina formatted fastq.gz file (sampleID_S01_L001_R1.fastq.gz, for example). Other formats will cause errors. Do not include a pathname since this is also the prefix that will be used for naming output files.
* ``-f,--fastqDir``: this is the directory in which the sample fastq.gz files can be found.
* ``-i,--intermediateDir``: this needs to be on a FIFO-capable server that is *not* the same as the folder containing the sample fastq.gz file. This protects the fastq files from accidental deletion or corruption if an error occurs during alignment.
* ``-p,--starParams``: this text file contains customizable parameters for STAR alignment (STAR is incredibly customizable, and hardcoding parameters into the code would make it much harder to optimize the alignment).
* ``-r,--refDir``: this directory needs to contain the reference GFF/GTF file, the reference fasta file, and the STAR indexes generated from those files.
* ``-o,--outputDir``: this is the directory into which the alignment output files will be placed.
* ``-a,--array``: DO NOT USE this argument when running the module individually. It allows the workflow wrapper to parallelize alignment with a slurm job array.
* ``-h,--help``: prints an informational message and exits.

### Stringtie Quantification of STAR Alignments
The script for this step is A2.stringtie-quant.sh, and uses the [StringTie]("https://ccb.jhu.edu/software/stringtie/", "Title") tool from Johns Hopkins. For biocore use, the Mamba environment ``/data/biocore/programs/mamba-envs/biocore-rna`` is hardcoded into the script and contains all the necessary programs.

As with the previous step, running this script independently can only be done with one sample per call, to allow parallelized efficiency.

#### Usage:
    sbatch A2.stringtie-quant.sh \
      --sampleID "$i" \
      --alignmentDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment" \
      --quantDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-stringtie" \
      --refGTF "/data/gencore/databases/reference_genomes/massilia/Massilia.MET4.gff"

#### Arguments:
* ``-s,--sampleID``: sample ID includes all text before the first underscore in a standard Illumina formatted fastq.gz file (sampleID_S01_L001_R1.fastq.gz, for example). Other formats will cause errors.
* ``-a,--alignmentDir``: the directory in which the BAM file for the specified sampleID is located.
* ``-q,--quantDir``: the directory into which the stringtie output will be placed. A subfolder will be created named for the provided sampleID.
* ``-r,--refGTF``: full pathname to the reference GTF or GFF file used for the STAR alignment.
* ``-l,--list``: DO NOT USE this argument when running the module individually. It allows the workflow wrapper to parallelize quantification with a slurm job array.
* ``-h,--help``: prints an informational message and exits.

### Generate Gene Count Matrix
The script(s) (A3.merge-quant.sh and 4 auxiliary scripts) for this step take the individual stringtie quant files as input and merge them to generate a complete gene count matrix for the dataset. The stringtie subfolders (named with each sample ID) should all be present in the same directory. In addition, summary QC statistics for all STAR alignments are created, and those alignments should all be in the same directory as each other. This script uses some auxiliary scripts, which are included in the source directory for this module. These are A4.prepDE.py (original source [here]("https://github.com/iandriver/RNA-sequence-tools/tree/master","Title")), A5.stringtie_expression_matrix.pl (original source [here]("https://github.com/griffithlab/rnaseq_tutorial/tree/master/scripts","Title")), A6.clean-counts.py (written in-house), and A7.star-summary.py (also written in-house).

#### Usage:
    sbatch A3.merge-quant.sh \
      --alignmentDir /path/to/alignment/output/directory \
      --quantDir /path/to/stringtie/output/directory \
      --scriptDir /path/to/auxiliary-script/directory \

#### Arguments:
* ``-s,--scriptDir``: the directory in which the four auxiliary scripts can be found.
* ``-q,--quantDir``: the directory in which the stringtie subfolders can be found.
* ``-a,--alignmentDir``: the directory in which the STAR alignment final logs can be found, for QC stats.
* ``-h,--help``: prints an informational message and exits.

## Statistical Analysis with R

### DE R scripts
Here, a single wrapper script allows each of the provided R scripts to be run from the same standard command line call. Currently, the DE packages that can be run are DEseq2, edgeR, and NOIseq (hopefully, more will be added in the future!) In addition to the input files needed for each script, the script name needs to be submitted to the wrapper.

Several of these packages are not currently publicly available on Sol, so they will need to be installed in each user's personal R library.

Depending on the setup of the STAR reference GTF, you may need to use the transcript count matrix instead of the gene count matrix. It may be good to double check the matrix files to see which one has the geneID (or transcriptID) that matches the ID names in the GTF file, particularly if you're working with a new species.

R packages needed for all scripts:
* `dplyr`
* `optparse`
* `viridis`

An example of how to call this script from the wrapper:

    sbatch ./B1.DEG.rscripts.sh \
      --script "deseq2 edgeR NOIseq" \
      --directory /path/to/output/directory \
      --geneMatrix /path/to/gene.count.matrix.csv \
      --comparisons /path/to/comparisons.csv

Input parameters:
* `-s,--script`: the name(s) of the R script that should be called. Acceptable options include deseq2, edgeR, and NOIseq, and the script name is case-insensitive. You can submit more than one name by providing a quoted, space-delimited list of names.
* `-d,--directory`: the name of the directory to which output files should be written. If an absolute pathname is not provided for the gene count matrix or comparisons file, the script will look for them in this folder.
* `-g,--geneMatrix`: the name of the file (absolute or within `--directory`) that contains the gene count matrix.
* `-c,--comparisons`: the name of the file (absolute or within `--directory`) that contains the pairwise comparison data for differential analysis. A sample comparison file is included here; essentially, sample names should be row names and comparison names should be the column headers. Within each comparison, each sample in the control group is represented by a 1 and each sample in the experimental group is represented by a -1. Any unused samples for a given comparison are represented as 0. The same file format can be used for all three differential expression scripts.

### B2.deseq2.R

This script uses the R package [DEseq2](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-014-0550-8,"Title") to normalize counts and determine differentially expressed genes using logistic regression models to calculate expression change and the Wald test to establish significance. If you choose to run only a single differential expression tool, this is the best one for use cases where each group has at least two samples.

Additional R packages needed:
* `DESeq2`
* `tibble`
* `ggplot2`
* `RColorBrewer`
* `ggrepel`
* `pheatmap`
* `stringr`


### B3.edgeR.R

This script uses the R package [edgeR](https://academic.oup.com/bioinformatics/article/26/1/139/182458?login=true,"Title"), designed for low-replicate datasets, to normalize counts and determine differentially expressed genes. EdgeR models gene expression with a negative binomial distribution and assesses differential expression using an adapted Fisher's exact test. This script has the highest memory demand, so if time and memory are an issue it can be skipped.

Additional R packages needed:
* `edgeR`
* `pheatmap`

### B4.noiseq.R

This script uses the R package [NOISeq](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4666377/,"Title") to normalize counts and determine differentially expressed genes. NOIseq filters low-count features taking the experimental design into consideration and corrects for batch effect as part of normalization. For data sets with groups containing only a single sample, this is the best script to use for differential analysis as it was designed for that use case. However, in a recent comparison of DE R tools, this one was not among the top scorers, and I'll hopefully be replacing it soon.

Additional R packages needed:
* `NOISeq`
* `tidyverse`

### Merge DEG Lists

After DE gene lists have been generated, they should be merged with equivalent fields. The merge_de python script (written in-house) allows results from the three scripts to be merged based on the Log2 fold change value and the p-adjusted value. I typically run the script on Sol with the mamba environment ``/data/biocore/programs/mamba-envs/biocore-rna`` activated, since it has the necessary Python version and packages installed.

#### Usage:

``python /path/to/merge_de_both.py comparisons.csv log2FCcutoff``

#### Arguments

This script needs to be run in the folder containing the output DEG lists from all three RNAseq R packages from module B (in the next release, I'll update this script to include any/all tools that are used for analysis, instead of hard-coding in these specific three packages, and to allow different directories to be chosen.)

The arguments that need to be provided are the comparisons.csv file used for module B, and the comparison names in the header row need to match the names of the files (which should be automatic if you run the module B scripts as written), as well as the desired log2 fold change cutoff for inclusion (use 0 to include all genes)

#### Output

For each comparison, the script will output a merged csv file containing the p-adjusted and log2 fold change values from each tool used, and the averages of both p-adjusted and log2 fold change values for all tools where they were generated.

In addition, Venn diagrams showing the correspondence in detected differential genes between the three packages are created for each comparison.

### Annotate Gene Lists

Script C2 annotates gene lists by joining a gene information file containing a list of gene IDs in the first column with any other file containing a subset of those gene IDs in its first column. Currently the script can take either a list of files or a full directory as input, but it can't distinguish between file types and image files will cause problems.

This step is really only necessary if your gene info file contains columns of information that weren't included in the original count matrix csv file. For example, if you run analysis with only the IDs, you'll want to re-annotate the results with the gene names; similarly, if you had a file containing descriptions and functional terms for each gene, you could annotate with that.

#### Usage:

    bash /data/to/C2.annotate-DEGs.sh \
      --infoFile /path/to/geneInfo.tab \
      --inputDirectory /path/to/deglists \
      (--list "file1 file2 file3")

#### Output

An annotated version of each input DEG list file.
