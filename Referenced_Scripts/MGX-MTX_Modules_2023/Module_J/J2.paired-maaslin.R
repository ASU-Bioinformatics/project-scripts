install.packages('devtools')
library(devtools)
install_github('biobakery/MTX_model')

library(data.table)
library(MTXmodel)

setwd('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann90-paired_pipeline/pathways/')

meta <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/metadata.txt', data.table=F)
dna <- fread("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann90_pipeline/pathways/pathway-abundance/normalized_pathabundance.tsv", data.table=F)
rna <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann90-paired_pipeline/pathways/pathway-abundance/normalized_pathabundance.tsv', data.table=F)

rownames(dna) <- dna[,1]
dna <- as.matrix(dna[,-1])
dna <- dna[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(dna)),]
dna <- dna[,colSums(dna) > 0]
colnames(dna) <- sub('_nohost_SQP_L001_RC_001_Abundance-CPM','',colnames(dna))
dna <- dna[grep('\\|',rownames(dna)),]

rownames(rna) <- rna[,1]
rna <- as.matrix(rna[,-1])
rna <- rna[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(rna)),]
rna <- rna[,colSums(rna) > 0]
colnames(rna) <- sub('_nohost_SQP_L001_RC_001_Abundance-CPM','',colnames(rna))
rna <- rna[grep('\\|',rownames(rna)),]

rownames(meta) <- meta[,1]

# use the same model for RNA that was used for the DNA Maaslin run, for consistency

fit_data <- MTXmodel(
  rna, meta, 'output', transform = "LOG",
  fixed_effects = c('Diagnosis'),
  random_effects = c('Age','Sex','Race'),
  reference = "Diagnosis,Healthy",
  normalization = "NONE",
  standardize = FALSE,
  input_dnadata = dna
)

