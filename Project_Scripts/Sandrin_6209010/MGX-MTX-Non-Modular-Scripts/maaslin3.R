# MaAsLin3 is the newest version of MaAsLin as of August 2024
# it looks like it takes compositionality into account better than MaAsLin2
# so I'm going to try it out here

if (!require('kableExtra', character.only = TRUE)) {
  install.packages('kableExtra')
}

# Install devtools if not present
if (!require('devtools', character.only = TRUE)) {
  install.packages('devtools')
}

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.19")

BiocManager::install('optparse')
BiocManager::install('biomaRt')

# Install MaAsLin 3
library("devtools")
install_github("biobakery/maaslin3")

for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr', 'kableExtra', 'optparse')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

#### READ IN ARGUMENTS ####
option_list <- list(
  make_option(c("-m", "--metadata"), type="character",
              default="/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/Publication-Preparation/metadata.txt",
              help="path to metadata file"),
  make_option(c("-i", "--input"), type="character",
              default="/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/4.humann-assembly-free-analysis/1.metaphlan-classifications/2.taxonomic-classifications/merged_metaphlan_table_species.txt",
              help="taxonomic file name and pathway"),
  make_option(c("-p", "--pathway"), type="character",
              default="/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/4.humann-assembly-free-analysis/2.DNA_UniRef50/1.pathway-abundance/normalized_pathabundance_unstratified.tsv"),
  make_option(c("-k", "--kegg"), type="character",
              default="/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/8.DNA-abundance-profiling/ko_map_module_completeness.tab"),
  make_option(c("-o", "--output"), type="character",
              default="/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/Publication-Preparation/",
              help="directory to output MaAsLin3 results")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

#### DEFINE CUSTOM VARIABLES ####

metaFile <- opts$metadata
outputDir <- opts$output
metaphlanTaxonomy <- opts$input
unirefPathways <- opts$pathway
keggModules <- opts$kegg
#### SET DIRECTORY AND UPLOAD METADATA ####

setwd(outputDir)

metadata <- data.frame(read.delim(metaFile, sep="\t"))

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]

#### ANALYZE MAASLIN ####

# input and format pathway input
pa.input <- data.frame(read.delim(unirefPathways, sep="\t"))

rownames(pa.input) <- pa.input[,1]
pa.input <- pa.input[,-1]

pa.input <- as.data.frame(t(pa.input))
rownames(pa.input) <- rownames(metadata)
# input and format taxonomy input (removed 's_' prefix from taxon names and shortened sample names in text editor)
tax.input <- data.frame(read.delim(metaphlanTaxonomy, sep="\t"))

rownames(tax.input) <- tax.input[,1]
tax.input <- tax.input[,-1]

tax.input <- as.data.frame(t(tax.input))
rownames(tax.input) <- rownames(metadata)

kegg.input <- data.frame(read.delim(keggModules, sep="\t"))
rownames(kegg.input) <- kegg.input[,2]
kegg.input <- kegg.input[,-1]
kegg.input <- kegg.input[,-1]
kegg.input <- kegg.input[,-1]
kegg.input <- as.data.frame(t(kegg.input))
rownames(kegg.input) <- rownames(metadata)

# factor the categorical variables in metadata
metadata$DiseaseState <- factor(metadata$DiseaseState)
metadata$Diagnosis <- factor(metadata$Diagnosis)
metadata$Sex <- factor(metadata$Sex)
metadata$Race <- factor(metadata$Race)
metadata$DiseaseGroup <- factor(metadata$DiseaseGroup)

# this provides taxonomic abundance comparisons
fit_out <- maaslin3(input_data = tax.input,
                    input_metadata = metadata,
                    output = 'test.maaslin.taxa.output.1',
                    formula = '~ DiseaseGroup + reads',
                    reference = 'DiseaseGroup,Healthy',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1,
                    coef_plot_vars = NULL,
                    heatmap_vars = "DiseaseGroup All Crohn's Disease,
                                      DiseaseGroup All Ulcerative Colitis,
                                      DiseaseGroup Other Disease"
                    )
                    

# this provides pathway abundance comparisons
fit_out.pa50 <- maaslin3(input_data = pa.input,
                        input_metadata = metadata,
                        output = 'u50.maaslin.pathway.output',
                        formula = '~ DiseaseGroup + reads',
                       reference = 'DiseaseGroup,Healthy',
                        normalization = 'TSS',
                        transform = 'LOG',
                        augment = TRUE,
                        standardize = TRUE,
                        max_significance = 0.1,
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        max_pngs = 100,
                        cores = 1,
                       coef_plot_vars = "DiseaseGroup All Crohn's Disease,
                                        DiseaseGroup All Ulcerative Colitis",
                       heatmap_vars = "DiseaseGroup All Crohn's Disease,
                                        DiseaseGroup All Ulcerative Colitis,
                                        DiseaseGroup Other Disease"
                       )

# this provides kegg module abundance comparisons
fit_out.kegg <- maaslin3(input_data = kegg.input,
                       input_metadata = metadata,
                       output = 'CAT.kegg.maaslin.output',
                       formula = '~ DiseaseGroup + reads',
                       reference = 'DiseaseGroup,Healthy',
                       normalization = 'TSS',
                       transform = 'LOG',
                       augment = TRUE,
                       standardize = TRUE,
                       max_significance = 0.1,
                       median_comparison_abundance = TRUE,
                       median_comparison_prevalence = FALSE,
                       max_pngs = 100,
                       cores = 1,
                       coef_plot_vars = "DiseaseGroup All Crohn's Disease,
                                        DiseaseGroup All Ulcerative Colitis",
                       heatmap_vars = "DiseaseGroup All Crohn's Disease,
                                        DiseaseGroup All Ulcerative Colitis,
                                        DiseaseGroup Other Disease")
