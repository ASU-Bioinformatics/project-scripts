# Module B
Module B contains the scripts needed to run the complete [KRAKEN2](https://github.com/DerrickWood/kraken2/wiki,"Title") pathway from taxonomic classification of reads through to ecological diversity calculations and visualization. The last script is in [R](https://www.r-project.org/,"Title") and is designed to be run manually to allow for customization of graphic parameters.

## B1.kraken.sh
The kraken2 installation is hardcoded into this script and refers to a copy on ``/data/biocore`` on Sol. If you have access to this drive, nothing else is needed. This script will use the Kraken2 software to classify reads based on k-mer similarity to sequences in the reference database. Input parameters include the directory containing the fastq.gz read files, the output directory where kraken2 results should be stored, and the directory containing the kraken2 directory the script should use.

A sample call:

    sbatch ./B1.kraken.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline \
      --dbDir /data/gencore/databases/kraken/k2_pluspf

The options are as follows; all are required:
* ``-i,--inputDir``: the directory containing the pre-processed fastq.gz files
* ``-o,--outDir``: the directory where kraken outputs will be stored
* ``-d,--dbDir``: the directory where the desired kraken database is stored

Output from this script includes raw kraken files, human readable kraken reports, fastq files containing classified reads, and fastq files containing unclassified reads. The script will automatically place all files in subfolders by output type.

## B2.bracken.sh
Since Kraken classifies reads to multiple taxonomic levels, relative abundance can't be easily determined from the output. [Bracken](https://github.com/jenniferlu717/Bracken,"Title") uses Bayesian methods to estimate abundance of each taxon at each taxonomic level in each sample, using the Kraken classification reports as input.

As with Kraken2, an installation of Bracken on ``/data/biocore`` is hardcoded into the script.

A sample call:

    sbatch ./B2.bracken.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/kraken-outputs/kraken-reports \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs \
      --dbDir /data/gencore/databases/kraken/k2_pluspf \
      --levels "S G F O C P K D"

The options are as follows; all are required except levels (the default taxonomic levels included in the analysis are shown in the call above):
* ``-i,--inputDir``: the directory containing the Kraken report files
* ``-o,--outDir``: the directory where Bracken results should be written
* ``-d,--dbDir``: the directory where the desired Kraken database is located
* ``-t,--levels``: a list containing a single letter for each taxonomic level for which analysis is desired (S=species, G=genus, F=family, O=order, C=class, P=phylum, K=kingdom, D=domain)

## B3.bracken-diversity.sh
This script uses the diversity tools from the helper software [KrakenTools](https://github.com/jenniferlu717/KrakenTools/tree/master/DiversityTools,"Title") to calculate alpha diversity for each sample and a beta diversity matrix for the population. An installation of KrakenTools on ``/data/biocore`` is hardcoded into the script and a conda environment (also on ``/data/biocore``) is used to provide the necessary Python packages.

A sample call:

    sbatch ./B3.bracken-diversity.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-raws \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-diversity

The options are as follows:
* ``-i,--inputDir``: the directory containing the raw Bracken output files
* ``-o,--outDir``: the directory where Bracken diversity results should be written

This script outputs a text file of Shannon, Simpson, Fisher, Berger-Parker, and Inverse Simpson alpha diversity metrics for each sample at each taxonomic level for which Bracken raw files were generated. For each taxonomic level, the values for each sample are condensed into summary files. A single text file containing the beta diversity matrix is also created.

## B4.krona-plots.sh
[Krona](https://github.com/marbl/Krona/wiki,"Title") creates interactive sunburst visualizations of taxonomic abundance that can be opened in an internet browser. It is a quick way to make useful and beautiful graphics for the complexity of multiple levels of taxonomic classification within a sample.
This script uses KrakenTools and [KronaTools](https://github.com/marbl/Krona/wiki/KronaTools,"Title") to convert Bracken reports into Krona plots with accompanying tabular data.

A sample call:

    sbatch ./B4.krona-plots.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-reports \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/krona-outputs \
      --htmlDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/krona-htmls \
      --krakenBin /data/biocore/programs/KrakenTools-1.2 \
      --kronaBin /data/biocore/programs/KronaTools-2.8.1/bin

The options are as follows. krakenBin and kronaBin are optional, and the values shown here are the defaults.
* ``-i,--inputDir``: the directory containing the Bracken report files
* ``-o,--outputDir``: the directory to which tabular Krona output should be written
* ``-h,--htmlDir``: the directory to which the html Krona plots should be written
* ``-k,--krakenBin``: the directory for the installation of KrakenTools
* ``-n,--kronaBin``: the bin directory within the installation KronaTools

The output from this script includes a folder of tabular classification and abundance data for each sample and a folder containing an HTML Krona plot for each sample.

## B5.diversity-visualization.R
This script takes the alpha and beta diversity metrics calculated in B3 and creates figures to visualize the data. It requires a metadata file to facilitate comparisons between sample groups.

This script is intended to be run manually, so graphics can be adjusted depending on sample size. It can be run for any taxonomic level, although the default in the script is for species level classifications.

Specialized R packages used include [ggplot2](https://ggplot2.tidyverse.org/,"Title"), [ggstatsplot](https://indrajeetpatil.github.io/ggstatsplot/,"Title"), and [corrplot](https://cran.r-project.org/web/packages/corrplot/vignettes/corrplot-intro.html,"Title").

Inputs needed include:
* tabular file containing all alpha diversity metrics and all metadata information for all samples. SampleID should be the name of the column containing the sample names, and the alpha diversity columns should be named Shannon, Berger.Parker, Simpson, Fisher, and InverseSimpson. The order of the columns isn't important, but metadata column names need to be added to the script manually.
    * The script expects this input file to be called ``species_alpha-diversity_forR.txt``
* tabular file containing the beta diversity matrix. B3 outputs a file containing additional information, so that information should be removed and any special characters cleaned from the sample names.
    * The script expects this input file to be called ``beta-diversity_forR.txt``

For alpha diversity, each categorical variable in the metadata is plotted as a box-violin plot for each diversity metric, and each quantitative variable is plotted as a scatter plot. The box-violin plots included statistical pairwise comparisons between groups, and will provide p-values for statistically significant differences. Linear regression can be used to evaluate the scatter plot if a linear relationship appears to be present.

For beta diversity, a correlation plot is used to visualize similarities between samples.

## Before Moving On

Check the Krona plots for every sample to make sure all the abundance estimates seem reasonable for your sample type (i.e., look for possible contamination). Any extreme outliers in terms of population abundance can be removed from downstream analysis if necessary, but don't be eager to eliminate samples.
