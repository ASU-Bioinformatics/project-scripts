#### code for normalizing metatranscriptome data by pathway and taxonomic classification ####
# input is humann/metaphlan tables
# can be found at https://github.com/zji90/microbiome_taxacontrol
# this file formats/prepares the data for the downstream analysis

library(data.table)
setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/assembly-free-taxonomy/humann3_uniref50")

#meta <- fread('data/metadata/hmp2_metadata.csv',data.table = F)
meta <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/metadata.txt', data.table=F)
#d <- fread('data/metagenome/pathabundances_3.tsv',data.table = F)
d <- fread("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/assembly-free-taxonomy/humann3_uniref50/pathways/pathway-abundance/normalized_pathabundance.tsv", data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
d <- d[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(d)),]
d <- d[,colSums(d) > 0]
colnames(d) <- sub('_SQP_L001_RC_001_Abundance-CPM','',colnames(d))
d <- d[grep('\\|',rownames(d)),]
gend <- log2(d+1)

rd <- sub('_CAG_[0-9]*$','',rownames(gend))
gend <- t(sapply(unique(rd),function(i) {
  colMeans(gend[rd==i,,drop=F])
}))

#### waiting for RNA metaphlan to complete ####

#d <- fread('data/metatranscriptome/pathabundances_3.tsv',data.table = F)
d <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/assembly-free-taxonomy/humann3_uniref50_paired/pathways/pathway-abundance/normalized_pathabundance.tsv', data.table=F)
rownames(d) <- d[,1]
d <- as.matrix(d[,-1])
d <- d[!grepl('^UNMAPPED|^UNINTEGRATED|unclassified',rownames(d)),]
d <- d[,colSums(d) > 0]
colnames(d) <- sub('_SQP_L001_RC_001_Abundance-CPM','',colnames(d))
d <- d[grep('\\|',rownames(d)),]
trad <- log2(d+1)

rd <- sub('_CAG_[0-9]*$','',rownames(trad))
trad <- t(sapply(unique(rd),function(i) {
  colMeans(trad[rd==i,,drop=F])
}))

rm('d')
int <- list(intersect(rownames(gend),rownames(trad)),intersect(colnames(gend),colnames(trad)))

gend <- gend[int[[1]],int[[2]]]
trad <- trad[int[[1]],int[[2]]]

#p <- fread('data/metagenome/taxonomic_profiles_3.tsv',data.table = F)
p <- fread('/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/assembly-free-taxonomy/humann3_uniref50/taxonomic-classifications/merged_metaphlan_table.txt', data.table=F)
rownames(p) <- p[,1]
p <- as.matrix(p[,-1])
colnames(p) <- sub('_SQP_L001_RC_001_metaphlan_bugs_list','',colnames(p))
p <- p[!grepl('^UNKNOWN',rownames(p)),]
p <- p[,colSums(p) > 0]

p <- p[,colnames(trad)]
p <- p[grep('k.*\\|p.*\\|c.*\\|o.*\\|f.*\\|g.*\\|s.*',rownames(p)),]
rd <- sub('_CAG_[0-9]*$','',rownames(p))
p <- t(sapply(unique(rd),function(i) {
  colSums(p[rd==i,,drop=F])
}))
p <- p/100

species <- sub('.*\\|','',rownames(trad))
species <- sub('\\.','|',species)
path <- sub('\\|.*','',rownames(trad))
rownames(p) <- sub('.*\\|g__','g__',rownames(p))

meta <- meta[match(colnames(gend),meta[,2]),]

res <- list(dna=gend,rna=trad,abu=p,species=species,pathway=path,meta=meta)
saveRDS(res,file='/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/assembly-free-taxonomy/humann3_uniref50_paired/proc.rds')

