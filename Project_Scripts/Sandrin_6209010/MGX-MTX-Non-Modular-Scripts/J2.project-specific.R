install.packages('devtools')
library(devtools)
install_github('biobakery/MTX_model')

library(data.table)
library(MTXmodel)

setwd('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/assembly-free-taxonomy/humann3_uniref90_paired/pathways/')

meta <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/metadata.txt', data.table=F)
dna <- fread("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/assembly-free-taxonomy/humann3_uniref90/pathways/pathway-abundance/normalized_pathabundance.tsv", data.table=F)
rna <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/assembly-free-taxonomy/humann3_uniref90_paired/pathways/pathway-abundance/normalized_pathabundance.tsv', data.table=F)

rownames(dna) <- dna[,1]
dna <- as.matrix(dna[,-1])
dna <- dna[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(dna)),]
dna <- dna[,colSums(dna) > 0]
colnames(dna) <- sub('_SQP_L001_RC_001_Abundance-CPM','',colnames(dna))
dna <- dna[grep('\\|',rownames(dna)),]

rownames(rna) <- rna[,1]
rna <- as.matrix(rna[,-1])
rna <- rna[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(rna)),]
rna <- rna[,colSums(rna) > 0]
colnames(rna) <- sub('_SQP_L001_RC_001_Abundance-CPM','',colnames(rna))
rna <- rna[grep('\\|',rownames(rna)),]

rownames(meta) <- meta[,1]

fit_data_1 <- MTXmodel(
  rna, meta, 'output', transform = "LOG",
  fixed_effects = c('Diagnosis', 'DiseaseState', 'Age'),
  reference = "Diagnosis,Healthy",
  normalization = "NONE",
  standardize = FALSE,
  input_dnadata = dna
)
