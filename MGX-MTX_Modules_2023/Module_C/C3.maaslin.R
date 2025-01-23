#### INSTALL AND LOAD MAASLIN ####

#if(!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")

#install.packages('remotes')
#remove.packages("Maaslin2")
#purge("Maaslin2")
#remotes::install_github(repo="biobakery/Maaslin2", force=TRUE)

#install.packages("corrr")
#install.packages("ggcorrplot")
#install.packages("FactoMineR")
#install.packages("factoextra")

pkg1 <- c("corrr", "ggcorrplot", "FactoMineR", "factoextra", "Maaslin2")
lapply(pkg1, library, character.only = TRUE)

#### DEFINE CUSTOM VARIABLES ####

# this is only for nonpaired data. Use MTX model scripts for paired metatranscriptomic data (section J)
# to prepare abundance file for R, remove hash from Pathway column and shorten the sample names

metaDir <- "/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/"
metaFile <- paste0(metaDir,"metadata.txt")
projectDir <- "/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann50_pipeline/"
outputDir <- paste0(projectDir, "pathways/differential-pathways")
inputDir <- paste0(projectDir, "pathways/pathway-abundance/")
normalizedAbundance <- paste0(inputDir, "normalized_pathabundance_unstratified.tsv")

#### SET DIRECTORY AND UPLOAD METADATA ####

setwd(outputDir)

metadata <- data.frame(read.delim(metaFile, sep="\t"))

rownames(metadata) <- metadata[,1]

#### ANALYZE MAASLIN ####

pa.input <- data.frame(read.delim(normalizedAbundance, sep="\t"))

rownames(pa.input) <- pa.input[,1]
pa.input <- pa.input[,-1]

#correct column names where needed (eg., if they start with a numeric character)
colnames(pa.input) <- c("MB-001", "MB-003", "MB-004", "MB-005", "MB-006", "MB-011", "MB-012", "MB-014", "MB-016", "MB-017",
                        "MB-018", "MB-020", "MB-021", "MB-023", "MB-024", "MB-025", "MB-028", "MB-033", "MB-035", "MB-037",
                        "MB-038", "MB-039", "MB-040", "MB-044", "MB-045", "MB-047", "MB-049", "MB-050", "MB-053", "MB-055")

# fixed effects are variables we think should have some sort of effect on the results
# we can also include random effects - such as sequencing batch, or an individual across a time course

fit_data4 = Maaslin2(
  input_data = pa.input,
  input_metadata = metadata,
  output = "maaslin2_90_model4",
  min_abundance = 0.0001,
  fixed_effects = c('Diagnosis'),
  random_effects = c('Age','Sex','Race'),
  reference = c("Diagnosis,Healthy")
  )

#### PCA ####

num <- length(rownames(pa.input))
input.mapped <- pa.input[3:num,] #this is to exclude unmapped and unintegrated reads
scaled.pa.input <- scale(input.mapped)
scaled.matrix <- cor(scaled.pa.input)

#pdf("Correlation_Plot.pdf")
#ggcorrplot(input.matrix)
#dev.off()

pca.scaled.input <- princomp(scaled.matrix)
summary(pca.scaled.input)
pca.scaled.input$loadings[, 1:2]

pdf("PCA_Dimensions.pdf")
fviz_eig(pca.scaled.input, addlabels = TRUE, title = "Pathway Abundance Component Contribution")
dev.off()

pdf("PCA_Plot.pdf", width=12, height=12)
fviz_pca_var(pca.scaled.input, col.var = "black", title = "Pathway Abundance PCA Plot")
dev.off()
