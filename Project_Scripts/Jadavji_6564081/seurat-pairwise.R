library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(Matrix)
library(rjson)
library(cowplot)
library(RColorBrewer)
library(grid)
library(readbitmap)
library(hdf5r)
library(data.table)
library(ggrepel)

#### read in all samples and create Seurat Objects ####
MFG.03.41 <- Load10X_Spatial(data.dir = "/data/gencore/sftp/n_jadavji/6564081_Spatial/spaceranger-output/Nafisa_MFG-03-41/outs",
                             slice = "MFG.03.41")
MFG.10.26 <- Load10X_Spatial(data.dir = "/data/gencore/sftp/n_jadavji/6564081_Spatial/spaceranger-output/Nafisa_MFG-10-26/outs",
                             slice = "MFG.10.26")
MFG.10.62 <- Load10X_Spatial(data.dir = "/data/gencore/sftp/n_jadavji/6564081_Spatial/spaceranger-output/Nafisa_MFG-10-62/outs",
                             slice = "MFG.10.62")
MFG.97.40 <- Load10X_Spatial(data.dir = "/data/gencore/sftp/n_jadavji/6564081_Spatial/spaceranger-output/Nafisa_MFG-97-40/outs",
                             slice = "MFG.97.40")

MFG.03.41@project.name <- "MFG.03.41"
MFG.03.41$orig.ident <- factor("MFG.03.41")
MFG.03.41@meta.data$sex <- factor("male")
MFG.03.41@meta.data$status <- factor("control")

MFG.10.26@project.name <- "MFG.10.26"
MFG.10.26$orig.ident <- factor("MFG.10.26")
MFG.10.26@meta.data$sex <- factor("female")
MFG.10.26@meta.data$status <- factor("control")

MFG.10.62@project.name <- "MFG.10.62"
MFG.10.62$orig.ident <- factor("MFG.10.62")
MFG.10.62@meta.data$sex <- factor("male")
MFG.10.62@meta.data$status <- factor("dementia")

MFG.97.40@project.name <- "MFG.97.40"
MFG.97.40$orig.ident <- factor("MFG.97.40")
MFG.97.40@meta.data$sex <- factor("female")
MFG.97.40@meta.data$status <- factor("dementia")

cortex <- merge(MFG.03.41,
                c(MFG.10.26, MFG.10.62, MFG.97.40),
                add.cell.ids=c("0341","1026", "1062", "9740"))
cortex <- SetIdent(cortex, value = cortex$orig.ident)

cortex <- PercentageFeatureSet(cortex, "^MT-", col.name = "percent_mito")
cortex.filtered <- cortex[, cortex$nFeature_Spatial > 10 & cortex$percent_mito < 50]

cortex.normalized <- SCTransform(cortex.filtered, assay = "Spatial", verbose = TRUE, method = "poisson")

cortex.normalized <- RunPCA(cortex.normalized, assay = "SCT", verbose = FALSE)
cortex.normalized <- FindNeighbors(cortex.normalized, reduction = "pca", dims = 1:30)
cortex.normalized <- FindClusters(cortex.normalized, verbose = FALSE)
cortex.normalized <- RunUMAP(cortex.normalized, reduction = "pca", dims = 1:30)
cortex.normalized <- RunTSNE(cortex.normalized, reduction = "pca", dims = 1:30)

#### pairwise by sex ####
Idents(cortex.normalized) <- cortex.normalized@meta.data$sex

pdf("cortex.normalized.sex.umap.pdf", width=8, height=6)
DimPlot(cortex.normalized, reduction = "umap", group.by = "sex") +
  plot_annotation("UMAP Visualization by Sex")
dev.off()

pdf("cortex.normalized.sex.tsne.pdf", width=8, height=6)
DimPlot(cortex.normalized, reduction = "tsne", group.by = "sex") +
  plot_annotation("tSNE Visualization by Sex")
dev.off()

#### pairwise by status ####
Idents(cortex.normalized) <- cortex.normalized@meta.data$status

pdf("cortex.normalized.status.umap.pdf", width=8, height=6)
DimPlot(cortex.normalized, reduction = "umap", group.by = "status") +
  plot_annotation("UMAP Visualization by Status")
dev.off()

pdf("cortex.normalized.status.tsne.pdf", width=8, height=6)
DimPlot(cortex.normalized, reduction = "tsne", group.by = "status") +
  plot_annotation("tSNE Visualization by Status")
dev.off()

#### integration ####
# create a list of the original data that we loaded to start with
st.list = list(MFG.03.41 = MFG.03.41, MFG.10.26 = MFG.10.26,
               MFG.10.62 = MFG.10.62, MFG.97.40 = MFG.97.40)

# purge any cells that fell to 0 counts at any point
st.list$MFG.03.41 <- st.list$MFG.03.41[, st.list$MFG.03.41$nCount_Spatial > 0]
st.list$MFG.10.26 <- st.list$MFG.10.26[, st.list$MFG.10.26$nCount_Spatial > 0]
st.list$MFG.10.62 <- st.list$MFG.10.62[, st.list$MFG.10.62$nCount_Spatial > 0]
st.list$MFG.97.40 <- st.list$MFG.97.40[, st.list$MFG.97.40$nCount_Spatial > 0]

# run SCT on all datasets
st.list = lapply(st.list, SCTransform, assay = "Spatial")

# need to set maxSize for PrepSCTIntegration to work
options(future.globals.maxSize = 2000 * 1024^2)  # set allowed size to 2K MiB

int.anchors <- FindIntegrationAnchors(st.list)
cortex.integrated <- IntegrateData(anchorset = int.anchors, dims = 1:20)

DefaultAssay(cortex.integrated) <- "integrated"

# Run the standard workflow for visualization and clustering
cortex.integrated <- ScaleData(cortex.integrated, verbose = FALSE)
cortex.integrated <- RunPCA(cortex.integrated, npcs = 30, verbose = FALSE)
# UMAP and Clustering
cortex.integrated <- RunUMAP(cortex.integrated, reduction = "pca", dims = 1:20)
cortex.integrated <- FindNeighbors(cortex.integrated, reduction = "pca", dims = 1:20)
cortex.integrated <- FindClusters(cortex.integrated, resolution = 0.5)

p1 <- DimPlot(cortex.integrated, reduction = "umap", group.by = "sex")
p2 <- DimPlot(cortex.integrated, reduction = "umap", group.by = "status")
pdf("cortex.integrated.umap.sex-and-status.pdf", width=8, height=4)
DimPlot(cortex.integrated, reduction = "umap", group.by = c("sex", "status")) +
  plot_annotation("UMAP Visualization of Integrated Data by Sex (left) and Status (right)")
dev.off()

#### differential markers ####

cortex.integrated <- SetIdent(cortex.integrated, value = "sex")
cortex.integrated <- PrepSCTFindMarkers(cortex.integrated, assay="SCT")
cortex.integrated.sex.markers <- FindAllMarkers(cortex.integrated, assay="SCT", )

cortex.integrated <- SetIdent(cortex.integrated, value = "status")
cortex.integrated <- PrepSCTFindMarkers(cortex.integrated, assay="SCT")
cortex.integrated.status.markers <- FindAllMarkers(cortex.integrated, assay="SCT", )

#### MA plots by sex ####
options(ggrepel.max.overlaps = 20)
cortex.integrated.sex.markers$log10padj <- -log10(cortex.integrated.sex.markers$p_val_adj)
cortex.integrated.sex.markers$avg_expr <- (cortex.integrated.sex.markers$pct.1 + cortex.integrated.sex.markers$pct.2)/2

for (i in c("female", "male")) {
  pdf(paste0("cortex.integrated.maplot.sex.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.sex.markers[cortex.integrated.sex.markers$cluster==i, ],
               mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=subset(cortex.integrated.sex.markers[cortex.integrated.sex.markers$cluster==i, ],
                                      (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                          aes(avg_expr, avg_log2FC, label=gene)) +
          ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Included") +
          xlab("Average Percent of Spots Expressing Gene") +
          ylab(paste0("Average Fold Change Between ", i, " and Others")))
  dev.off()
}

cortex.integrated.sex.markers.nomito <- cortex.integrated.sex.markers[!grepl('^MT-', cortex.integrated.sex.markers$gene),]

for (i in c("MFG.03.41", "MFG.10.26", "MFG.10.62", "MFG.97.40")) {
  pdf(paste0("cortex.integrated.maplot.sex.nomito.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.sex.markers.nomito[cortex.integrated.sex.markers.nomito$cluster==i, ],
               mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=subset(cortex.integrated.sex.markers.nomito[cortex.integrated.sex.markers.nomito$cluster==i, ],
                                      (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                          aes(avg_expr, avg_log2FC, label=gene)) +
          ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Excluded") +
          xlab("Average Percent of Spots Expressing Gene") +
          ylab(paste0("Average Fold Change Between ", i, " and Others")))
  dev.off()
}

