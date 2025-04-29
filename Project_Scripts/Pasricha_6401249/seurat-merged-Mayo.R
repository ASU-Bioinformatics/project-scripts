library(dplyr)
library(Seurat)
library(patchwork)

#### read in all samples and create Seurat Objects ####
setwd("/data/gencore/analysis_projects/6401249_Pasricha_scRNAseq/globus-data/092022-2A/outs")

A.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
A <- CreateSeuratObject(counts = A.data, project = "092022.2A", min.cells = 3, min.features = 200)

setwd("/data/gencore/analysis_projects/6401249_Pasricha_scRNAseq/globus-data/092022-2B/outs")

B.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
B <- CreateSeuratObject(counts = B.data, project = "092022.2B", min.cells = 3, min.features = 200)

setwd("/data/gencore/analysis_projects/6401249_Pasricha_scRNAseq/globus-data/092022-2C/outs")

C.data <- Read10X(data.dir = "filtered_feature_bc_matrix")
C <- CreateSeuratObject(counts = C.data, project = "092022.2C", min.cells = 3, min.features = 200)

#move to output directory
setwd("/data/gencore/analysis_projects/6401249_Pasricha_scRNAseq/seurat-merged-output")

#### merge Seurat Objects ####

M <- merge(A, y = c(B, C),
           add.cell.ids = c("092022.2A", "092022.2B", "092022.2C"),
           project = "PasrichaAllSamples")

#### find and set cutoff thresholds using plots ####
M[["percent.mt"]] <- PercentageFeatureSet(M, pattern = "^MT-")
#M <- subset(M, subset = nFeature_RNA > 200 & nFeature_RNA < 20000 & percent.mt < 20)

#normalize the data
M <- NormalizeData(M)

#determine percent mitochondrial
#pdf("qcmetrics-M.pdf", height=7, width=14)
#VlnPlot(M, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

#subset by quality
M.lowMT <- subset(M, subset = nFeature_RNA > 200 & nFeature_RNA < 9000 & percent.mt < 10)
#pdf("qcmetrics-M.lowMT.pdf", height=7, width=14)
#VlnPlot(M.lowMT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
#dev.off()

#normalize the data
M.lowMT <- NormalizeData(M.lowMT)

#identify highly variable features
M.lowMT <- FindVariableFeatures(M.lowMT, selection.method = "vst", nfeatures = 4000)
top10 <- head(VariableFeatures(M.lowMT), 10)
plot1 <- VariableFeaturePlot(M.lowMT)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
#pdf("variablefeatures-top10-M.lowMT.pdf", width=20, height=7)
#plot1 + plot2
#dev.off()

#scale data
all.genes <- rownames(M.lowMT)
M.lowMT <- ScaleData(M.lowMT, features = all.genes)

#linear dimensional reduction (PCA)
M.lowMT <- RunPCA(M.lowMT, features = VariableFeatures(object = M.lowMT))
VizDimLoadings(M.lowMT, dims = 1:4, reduction = "pca")
pdf("PCA-by-samples.M.lowMT.pdf")
DimPlot(M.lowMT, reduction = "pca")
dev.off()

pdf("PC.6dim.heatmaps.M.lowMT.pdf")
heatmap <- DimHeatmap(M.lowMT, dims = c(1,2,3,4,5,6), cells = 500, balanced = TRUE, combine=TRUE)
dev.off()

# determine dimensionality of dataset
ElbowPlot(M.lowMT, ndims=50)

# this dataset doesn't have a clear elbow, more like a gentle curve. Try running with 33 PC's
# 42 PC's resulted in a lot of singletons, while 13 PC's resulted in almost no visible separation between cell clusters
# higher resolutions led to cell clusters that were highly overlapping

#### run with calculated number of PCs ####
M.lowMT.PC33 <- FindNeighbors(M.lowMT, dims = 1:33)
M.lowMT.PC33 <- FindClusters(M.lowMT.PC33, resolution = 0.4)

# run UMAP/tSNE
M.lowMT.PC33 <- RunUMAP(M.lowMT.PC33, dims = 1:33)
M.lowMT.PC33 <- RunTSNE(M.lowMT.PC33, dims = 1:33)
pdf("UMAP.M.lowMT.PC33.res0.4.pdf")
DimPlot(M.lowMT.PC33, reduction = "umap", label = TRUE)
dev.off()

M.lowMT.PC33 <- SetIdent(M.lowMT.PC33, value = "seurat_clusters")
pdf("TSNE.M.lowMT.PC33.res1.2.pdf")
DimPlot(M.lowMT.PC33, reduction = "tsne", label = TRUE)
dev.off()

saveRDS(M.lowMT.PC33, file = "M.lowMT.PC33.res0.4.rds")

# identify cluster biomarkers
M.lowMT.PC33.markers <- FindAllMarkers(M.lowMT.PC33, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
top2markers.M.lowMT.PC33.clusters <- M.lowMT.PC33.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

pdf("top2markers.lowMT.PC33.clusters0-5.pdf", width=14, height=7)
FeaturePlot(M.lowMT.PC33, features = top2markers.M.lowMT.PC33.clusters$gene[2:11])
dev.off()

pdf("top2markers.lowMT.PC33.clusters6-11.pdf", width=14, height=7)
FeaturePlot(M.lowMT.PC33, features = top2markers.M.lowMT.PC33.clusters$gene[12:23])
dev.off()

pdf("top2markers.M.lowMT.PC33.heatmap.pdf", width=28, height=14)
M.lowMT.PC33.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(M.lowMT.PC33, features = top10$gene) + NoLegend()
dev.off()



#### some more visuals ####

features <- unique(top2markers.M.lowMT.PC33.clusters$gene)

pdf("celltype.dotplot.lowMT.PC33.Res0.4.pdf", height=6, width=16)
DotPlot(M.lowMT.PC33, features = features) + RotatedAxis()
dev.off()

pdf("celltype.features.heatmap.lowMT.PC33.Res0.4.pdf", height=8, width=16)
DoHeatmap(subset(M.lowMT.PC33, downsample = 100), features = features, size = 3)
dev.off()

#### cell type proportion analysis between samples ####

#M.lowMT.PC33$CellType <- Idents(M.lowMT.PC33)

#devtools::install_github("rpolicastro/scProportionTest")

library("scProportionTest")
library("ggplot2")
M.utils <- sc_utils(M.lowMT.PC33)

# these are output as sample_2_vs_sample_1, where sample_1 is the base/control
# negative fold change means fewer cells in sample_2 than in sample_1
M.utils.BvA <- permutation_test(M.utils, cluster_identity = "seurat_clusters",
                                           sample_1 = "092022.2A", sample_2 = "092022.2B",
                                           sample_identity = "orig.ident")

M.utils.CvA <- permutation_test(M.utils, cluster_identity = "seurat_clusters",
                                           sample_1 = "092022.2A", sample_2 = "092022.2C",
                                           sample_identity = "orig.ident")

M.utils.CvB <- permutation_test(M.utils, cluster_identity = "seurat_clusters",
                                           sample_1 = "092022.2B", sample_2 = "092022.2C",
                                           sample_identity = "orig.ident")


plot.BvA <- permutation_plot(M.utils.BvA) +
  ggtitle("Differential Clustering in 092022.2B vs. 092022.2A")

plot.CvA <- permutation_plot(M.utils.CvA) +
  ggtitle("Differential Clustering in 092022.2C vs. 092022.2A")

plot.CvB <- permutation_plot(M.utils.CvB) +
  ggtitle("Differential Clustering in 092022.2C vs. 092022.2B")

pdf("differential-clustering-B_v_A.pdf", height=7, width=7)
print(plot.BvA)
dev.off()

pdf("differential-clustering-C_v_A.pdf", height=7, width=7)
print(plot.CvA)
dev.off()

pdf("differential-clustering-C_v_B.pdf", height=7, width=7)
print(plot.CvB)
dev.off()


### multi-group comparisons with CDTS ####
install.packages("hablar")
library("tidyverse")
library("hablar")
library("broom")
library("patchwork")

source("~/CTDS_function.R")
alphaDivScores <- CTDS.score(M.lowMT.PC33, sample = "orig.ident", cell.type = "seurat_clusters")
# increased CTDS score (closer to 0, farther from -1) implies a more uniform distribution of cell types #

#### differential expression by cluster ####
M.1.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "1", ident.2 = NULL, only.pos = TRUE)
M.2.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "2", ident.2 = NULL, only.pos = TRUE)
M.3.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "3", ident.2 = NULL, only.pos = TRUE)
M.4.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "4", ident.2 = NULL, only.pos = TRUE)
M.5.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "5", ident.2 = NULL, only.pos = TRUE)
M.6.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "6", ident.2 = NULL, only.pos = TRUE)
M.7.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "7", ident.2 = NULL, only.pos = TRUE)
M.8.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "8", ident.2 = NULL, only.pos = TRUE)
M.9.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "9", ident.2 = NULL, only.pos = TRUE)
M.10.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "10", ident.2 = NULL, only.pos = TRUE)
M.11.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "11", ident.2 = NULL, only.pos = TRUE)
M.0.de.markers <- FindMarkers(M.lowMT.PC33, ident.1 = "0", ident.2 = NULL, only.pos = TRUE)

# save workspace here #

degenes <- list("0" = M.0.de.markers,
                "1" = M.1.de.markers,
                "2" = M.2.de.markers,
                "3" = M.3.de.markers,
                "4" = M.4.de.markers,
                "5" = M.5.de.markers,
                "6" = M.6.de.markers,
                "7" = M.7.de.markers,
                "8" = M.8.de.markers,
                "9" = M.9.de.markers,
                "10" = M.10.de.markers,
                "11" = M.11.de.markers)

write.table(M.0.de.markers, file="deg.cluster0_vs_others.txt", sep="\t")
write.table(M.1.de.markers, file="deg.cluster1_vs_others.txt", sep="\t")
write.table(M.2.de.markers, file="deg.cluster2_vs_others.txt", sep="\t")
write.table(M.3.de.markers, file="deg.cluster3_vs_others.txt", sep="\t")
write.table(M.4.de.markers, file="deg.cluster4_vs_others.txt", sep="\t")
write.table(M.5.de.markers, file="deg.cluster5_vs_others.txt", sep="\t")
write.table(M.6.de.markers, file="deg.cluster6_vs_others.txt", sep="\t")
write.table(M.7.de.markers, file="deg.cluster7_vs_others.txt", sep="\t")
write.table(M.8.de.markers, file="deg.cluster8_vs_others.txt", sep="\t")
write.table(M.9.de.markers, file="deg.cluster9_vs_others.txt", sep="\t")
write.table(M.10.de.markers, file="deg.cluster10_vs_others.txt", sep="\t")
write.table(M.11.de.markers, file="deg.cluster11_vs_others.txt", sep="\t")

library(biomaRt)
mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

allGenes <- lapply(degenes,
                   function(x) {
                     as.list(getBM(filters = "external_gene_name",
                                   attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                                   values = rownames(x),
                                   mart = mart))
                   })
allGenes <- lapply(allGenes,
                   function(x) {
                     data.frame(x)})

idmaps <- sapply(names(allGenes),
                 function (z) {
                   degenes[[z]]$external_gene_name <- rownames(degenes[[z]])
                   merge(x = degenes[[z]], y = allGenes[[z]], by="external_gene_name", all.x=TRUE)
                 },
                 simplify = FALSE,
                 USE.NAMES = TRUE)

sapply(names(idmaps),
       function (x) {
         write.table(idmaps[[x]], file=paste0("annotated.deg.",x,".txt"), sep="\t")
       })
rownames(idmaps) <- paste(idmaps$ensembl_gene_id, idmaps$external_gene_name)


degenesRel <- lapply(degenes,
                  function(x) {
                    sorted <- x[order(x$p_val_adj),]
                    lfcCut <- sorted[sorted$avg_log2FC > 0.58,]
                    qval <- lfcCut[lfcCut$p_val_adj < 0.05, ]
                    return(qval)})




pkg2 <- c('igraph', 'clusterProfiler', 'enrichplot', 'ggnewscale', 'org.Mm.eg.db')
lapply(pkg2,
       function(x) library(x, character.only = TRUE))

degenes.list <- sapply(names(degenesRel),
                       function(x) {
                         enrichGO(gene = row.names(degenesRel[[x]]),
                                  OrgDb = 'org.Hs.eg.db',
                                  ont = "ALL",
                                  keyType = 'SYMBOL',
                                  pAdjustMethod = "BH",
                                  pvalueCutoff = 0.01,
                                  qvalueCutoff = 0.01)
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE)

## save workspace here ##

lapply(names(degenes.list),
       function (y) {
         if (!is.null(degenes.list[[y]])) {
           if (dim(degenes.list[[y]]@result)[[1]] != 0) {
             plot.bar <- mutate(degenes.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore")
             pdf(file = paste0("de-functions.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(degenes.list[[y]], showCategory=20)
             pdf(file = paste0("de-functions.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)

degenes.list <- sapply(names(degenes.list),
                         function (y) {
                           if (!is.null(degenes.list[[y]])) {
                             if (dim(degenes.list[[y]]@result)[[1]] > 0) {
                               pairwise_termsim(degenes.list[[y]])
                             }
                           }
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

lapply(names(degenes.list),
       function (y) {
         if (!is.null(degenes.list[[y]])) {
           if (dim(degenes.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(degenes.list[[y]],
                                   cex_label_category=0.7,
                                   cex_category=0.5,
                                   cex_line=0.5)
             pdf(file = paste0("de-functions.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

rd.degenes.list <- sapply(names(degenes.list),
                            function (x) {
                              if (!is.null(degenes.list[[x]])) {
                                setReadable(degenes.list[[x]],
                                            OrgDb = 'org.Hs.eg.db',
                                            keyType = "ENSEMBL")
                              }},
                            simplify = FALSE,
                            USE.NAMES = TRUE)

lapply(names(rd.degenes.list),
       function (x) {
         if (!is.null(rd.degenes.list[[x]])) {
           plot.cnet.gene <- cnetplot(rd.degenes.list[[x]],
                                      M = degenes.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(rd.degenes.list[[x]],
                                     M = degenes.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(rd.degenes.list[[x]],
                                     M = degenes.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("de-functions.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

sapply(names(degenes.list),
       function (x) {
         if (!is.null(rd.degenes.list[[x]])) {
           write.table(degenes.list[[x]]@result, file = paste0(x, ".clusterProfiler.sig.txt"), quote=FALSE)
         }})


