# Module C
The scripts in this model use [HUMAnN](https://huttenhower.sph.harvard.edu/humann/,"Title") to characterize the functional profile of an MGX or MTX data set using similarity to protein clusters from the [UniRef](https://www.uniprot.org/help/uniref,"Title") database. Along with this, [MetaPhlAn](https://huttenhower.sph.harvard.edu/metaphlan/,"Title") classifies reads into taxa so that each functional annotation can be paired with a taxonomic classification.

The scripts in this module use the Conda environment `/data/biocore/programs/conda-envs/humann-env/`, installed on the /data/biocore/ mount on Sol; if you need to create your own environment, you can use the provided yml files and replace the hardcoded pathnames for the environment in the scripts.

## C1.humann-metaphlan.wrapper.sh
This script performs the bulk of computational work in this module, and is where HUMAnN is called. HUMAnN requires concatenated fastq files, and the script will concatenate them if they haven't been already. In addition, this script can run either a standalone HUMAnN pipeline (for MGX or unpaired MTX data), or a paired HUMAnN pipeline where the MTX data uses the taxonomic population abundances calculated for their DNA counterparts.

The script is hardcoded to use a conda environment at ``/data/biocore/programs/conda-envs/humann-env/``, which can be changed to your own local conda environment if running off Sol.

A sample call for MGX data:

    bash ./C1.humann-metaphlan.wrapper.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --outputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann90_pipeline \
      --dbType u90

And a sample call for paired MTX data:

    bash ./C1.humann-metaphlan.wrapper.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq/concatenated \
      --outputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann90-paired_pipeline \
      --dnaRef /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann90_pipeline \
      --concat --paired --dbType u90

All available arguments are as follows:
* ``-i,--inputDir``: the directory containing the fastq files, whether they are already concatenated or not
* ``-o,--outDir``: the directory to which humann output should be written
* ``-r,--dnaRef``: the directory in which the DNA humann output is found for paired MTX analysis
* ``-d,--dbType``: this can be either u90 or u50. I have the database pathnames hardcoded into the script, and this serves as a toggle between them.
* ``-c,--concat``: if this flag is provided, the script looks for already concatenated fastq files
* ``-p,--paired``: if this flag is provided, the script runs a paired MTX analysis
* ``-l,--list``: if you don't want to run all the samples in the fastq directory, include a list of sample IDs following this tag. Sample IDs should include all text before "_SQP_L00x_Rx_001.fastq.gz" as this is the suffix the script looks for.

## C2.humann-manipulations.sh
This script does some helpful data cleanup and organization of the HUMAnN output. First, it normalizes all the pathway abundance and gene family tables into counts per million files, then merges those files to create a normalized matrix. Second, it reformats that data into a stratified and unstratified version. The stratified table shows the abundance of each pathway or gene family split by the taxonomic classifications associated with it (i.e., pathway 1 would be divided into pathway 1 from species A and pathway 1 from species B), while the unstratified table doesn't include the taxonomic information.

In addition to the data formatting, this script organizes the files into a structure convenient for data return.
This script is faster to run interactively than to set up as a script, so I typically run it manually. All that is needed to run it is the HUMAnN conda environment and the output directory from C1.

## C3.maaslin.R
[MaAsLin2](https://github.com/biobakery/Maaslin2,"Title") is an R package designed to model microbial populations in terms of multiple metadata groupings. Since this is an R script I run it manually for customization - especially customization of the model. Typically I will include any variable that I expect to influence variation, or that preliminary data suggests influences variation, as one of the fixed effects. The random effects are those factors that come from random sampling of a population, or where the span rather than the precise values could influence the model. For example, I've used factors such as Age, Race, and Sex as random effects, which smooths out  biological variation due to random sampling rather than to the factors of interest (such as Diagnosis or Treatment).

The model parameters used to run this script should be copied into a text file and saved in the project folder; the same parameters should be used for both UniRef50 and UniRef90 analysis.

In addition to the MaAsLin2 modeling, this script calculates principal components and generates a PCA plot based on pathway abundance.

Input for this script is a metadata text file and the normalized path abundance text file from HUMAnN. The normalized path abundance file does need to be prepared slightly, by removing the hash sign from the header row and renaming the sample names to their short IDs to match the metadata file.

## C4.metaphlan-heatmap.R
Using the file merged_metaphlan_table_species.txt created by the data formatting in C2, this script creates a clustered heatmap of metaphlan-predicted species across all samples. This provides a good visual for between-group differences in population fractions, as well as for whether samples in the same metadata categories cluster together or not.

This is also R so I run it manually. Any special characters in the file names will need to be corrected after important the species table.

## C5.metaphlan-alpha-diversity.R
This script relies on the R package [vegan](https://github.com/vegandevs/vegan,"Title") to calculate principal components and generate a PCA plot based on MetaPhlAn species-level taxonomic classifications, again using merged_metaphlan_table_species.txt as input.

Additionally, this script calculates, compares, and plots alpha diversity metrics for all samples, using the methods Shannon, Simpson, and InverseSimpson.

Outfiles include PCA plots for each metadata category and a box-violin plot for each metadata category for the three alpha diversity methods above.

## Before moving on
I always enjoy comparing the pathway abundance PCA plots to the MetaPhlAn PCA plots, to see how microbial population and pathway abundance cluster differently between the samples in the data set. For example, samples from similar environments across a time scale are more likely to have similar populations and different pathway abundances, while populations from different environments are more likely to have different populations.
