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
cortex <- PercentageFeatureSet(cortex, "^RP", col.name = "percent_ribosomal")

pdf("cortex.general-stats-violin.pdf", width=12, height=12)
VlnPlot(cortex, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf("cortex.general-stats-spatial.pdf", width=24, height=24)
SpatialFeaturePlot(cortex, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"), ncol = 4, pt.size.factor = 3)
dev.off()

pdf("cortex.general-stats-spatial-1062-only.pdf", width=24, height=24)
SpatialFeaturePlot(subset(x = cortex, idents = "MFG.10.62"),
                   features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
                   ncol = 4, pt.size.factor = 12)
dev.off()

#### filtered ####
# to recreate the violin plot after downstream analysis has been run, reset default assay and active idents
Idents(cortex.filtered) <- cortex.filtered$orig.ident
DefaultAssay(cortex.filtered) <- "Spatial"

# then, to reset them to the downstream state again
Idents(cortex.filtered) <- cortex.filtered$seurat_clusters
DefaultAssay(cortex.filtered) <- "SCT"

cortex.filtered <- cortex[, cortex$nFeature_Spatial > 10 & cortex$percent_mito < 50]
pdf("cortex.filtered.general-stats-violin.pdf", width=12, height=9)
VlnPlot(cortex.filtered, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend() +
  plot_annotation(title="General Metrics for All Samples", subtitle="Filtered to >10 features, <50% mitochondrial")
dev.off()


pdf("cortex.filtered.1062-only.general-stats-spatial.pdf", width=24, height=24)
SpatialFeaturePlot(subset(x = cortex.filtered, idents = "MFG.10.62"),
                   features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
                   ncol = 4, pt.size.factor = 12) +
  plot_annotation(title="General Metrics for All Samples", subtitle="Filtered to >10 features, <50% mitochondrial")
dev.off()

pdf("cortex.filtered.general-stats-spatial.pdf", width=24, height=24)
SpatialFeaturePlot(cortex.filtered,
                   features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
                   ncol = 4, pt.size.factor = 3) +
  plot_annotation(title="General Metrics for All Samples", subtitle="Filtered to >10 features, <50% mitochondrial")
dev.off()

#### SCTransform to Normalize ####
cortex.normalized <- SCTransform(cortex.filtered, assay = "Spatial", verbose = TRUE, method = "poisson")

cortex.normalized <- RunPCA(cortex.normalized, assay = "SCT", verbose = FALSE)
cortex.normalized <- FindNeighbors(cortex.normalized, reduction = "pca", dims = 1:30)
cortex.normalized <- FindClusters(cortex.normalized, verbose = FALSE)
cortex.normalized <- RunUMAP(cortex.normalized, reduction = "pca", dims = 1:30)
cortex.normalized <- RunTSNE(cortex.normalized, reduction = "pca", dims = 1:30)

pdf("cortex.normalized.umap.pdf", width=12, height=6)
DimPlot(cortex.normalized, reduction = "umap", group.by = c("ident", "orig.ident")) +
  plot_annotation("UMAP Visualization by Clusters (left) and Samples (right)")
dev.off()

pdf("cortex.normalized.tsne.pdf", width=12, height=6)
DimPlot(cortex.normalized, reduction = "tsne", group.by = c("ident", "orig.ident")) +
  plot_annotation("tSNE Visualization by Clusters (left) and Samples (right)")
dev.off()

pdf("cortex.normalized.spatial-clusters.pdf", width=12, height=12)
SpatialDimPlot(cortex.normalized, ncol = 2, pt.size.factor = 3) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster")
dev.off()

pdf("cortex.normalized.MFG1062-only.spatial-clusters.pdf", width=12, height=12)
SpatialDimPlot(subset(x = cortex.normalized, subset=orig.ident == "MFG.10.62"), ncol = 2, pt.size.factor = 12) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster")
dev.off()

pdf("cortex.normalized.spatial-clusters.alpha.pdf", width=12, height=12)
SpatialDimPlot(cortex.normalized, ncol = 2, pt.size.factor = 3, alpha=c(0.1,1)) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster")
dev.off()

pdf("cortex.normalized.MFG1062-only.spatial-clusters.alpha.pdf", width=12, height=12)
SpatialDimPlot(subset(x = cortex.normalized, subset=orig.ident == "MFG.10.62"),
               ncol = 2, pt.size.factor = 12, alpha=c(0.1,1)) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster")
dev.off()

pdf("cortex.normalized.spatial-clusters.MFG-03-41.pdf")
SpatialDimPlot(cortex.normalized,
               cells.highlight = CellsByIdentities(cortex.normalized),
               images = "MFG.03.41",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-03-41",
                  subtitle = "Calculated with normalized data from the merged sample set.")
dev.off()

pdf("cortex.normalized.spatial-clusters.MFG-10-26.pdf")
SpatialDimPlot(cortex.normalized,
               cells.highlight = CellsByIdentities(cortex.normalized),
               images = "MFG.10.26",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-10-26",
                  subtitle = "Calculated with normalized data from the merged sample set.")
dev.off()

pdf("cortex.normalized.spatial-clusters.MFG-10-62.pdf")
SpatialDimPlot(cortex.normalized,
               cells.highlight = CellsByIdentities(cortex.normalized),
               images = "MFG.10.62",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 12) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-10-62",
                  subtitle = "Calculated with normalized data from the merged sample set.")
dev.off()

pdf("cortex.normalized.spatial-clusters.MFG-97-40.pdf")
SpatialDimPlot(cortex.normalized,
               cells.highlight = CellsByIdentities(cortex.normalized),
               images = "MFG.97.40",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-97-40",
                  subtitle = "Calculated with normalized data from the merged sample set.")
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

p1 <- DimPlot(cortex.integrated, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(cortex.integrated, reduction = "umap", label = TRUE)
pdf("cortex.integrated.umap.cluster-and-sample.pdf", width=8, height=4)
DimPlot(cortex.integrated, reduction = "umap", group.by = c("ident", "orig.ident")) +
  plot_annotation("UMAP Visualization of Integrated Data by Clusters (left) and Samples (right)")
dev.off()

pdf("cortex.integrated.umap.clusters-by-sample.pdf", width=16, height=4)
DimPlot(cortex.integrated, reduction = "umap", split.by = "orig.ident") +
  plot_annotation(title="Cluster Representation in Each Sample")
dev.off()

#### spatial clusters plots ####
pdf("cortex.integrated.spatial-clusters.pdf", width=12, height=12)
SpatialDimPlot(cortex.integrated, ncol = 2, pt.size.factor = 3) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster", subtitle = "Integration-corrected Normalized Data")
dev.off()

pdf("cortex.integrated.MFG1062-only.spatial-clusters.pdf", width=12, height=12)
SpatialDimPlot(subset(x = cortex.integrated, subset=orig.ident == "MFG.10.62"), ncol = 2, pt.size.factor = 12) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster", subtitle = "Integration-corrected Normalized Data")
dev.off()

pdf("cortex.integrated.spatial-clusters.alpha.pdf", width=12, height=12)
SpatialDimPlot(cortex.integrated, ncol = 2, pt.size.factor = 3, alpha=c(0.1,1)) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster", subtitle = "Integration-corrected Normalized Data")
dev.off()

pdf("cortex.integrated.MFG1062-only.spatial-clusters.alpha.pdf", width=12, height=12)
SpatialDimPlot(subset(x = cortex.integrated, subset=orig.ident == "MFG.10.62"),
               ncol = 2, pt.size.factor = 12, alpha=c(0.1,1)) +
  plot_annotation("Spatial Visualization of Each Sample by Cluster", subtitle = "Integration-corrected Normalized Data")
dev.off()

pdf("cortex.integrated.spatial-clusters.MFG-03-41.pdf")
SpatialDimPlot(cortex.integrated,
               cells.highlight = CellsByIdentities(cortex.integrated),
               images = "MFG.03.41",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-03-41",
                  subtitle = "Calculated with integrated data from the merged sample set.")
dev.off()

pdf("cortex.integrated.spatial-clusters.MFG-10-26.pdf")
SpatialDimPlot(cortex.integrated,
               cells.highlight = CellsByIdentities(cortex.integrated),
               images = "MFG.10.26",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-10-26",
                  subtitle = "Calculated with integrated data from the merged sample set.")
dev.off()

pdf("cortex.integrated.spatial-clusters.MFG-10-62.pdf")
SpatialDimPlot(cortex.integrated,
               cells.highlight = CellsByIdentities(cortex.integrated),
               images = "MFG.10.62",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 12) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-10-62",
                  subtitle = "Calculated with integrated data from the merged sample set.")
dev.off()

pdf("cortex.integrated.spatial-clusters.MFG-97-40.pdf")
SpatialDimPlot(cortex.integrated,
               cells.highlight = CellsByIdentities(cortex.integrated),
               images = "MFG.97.40",
               facet.highlight = TRUE,
               ncol = 5, pt.size.factor = 3) +
  plot_annotation(title = "Spatial Visualization of Each Cluster in MFG-97-40",
                  subtitle = "Calculated with integrated data from the merged sample set.")
dev.off()

#### differential markers ####

cortex.integrated <- SetIdent(cortex.integrated, value = "orig.ident")
cortex.integrated <- PrepSCTFindMarkers(cortex.integrated, assay="SCT")
cortex.integrated.sample.markers <- FindAllMarkers(cortex.integrated, assay="SCT", )

cortex.integrated <- SetIdent(cortex.integrated, value = "seurat_clusters")
cortex.integrated <- PrepSCTFindMarkers(cortex.integrated, assay="SCT")
cortex.integrated.cluster.markers <- FindAllMarkers(cortex.integrated, assay="SCT", )

#### MA plots including mitochondrial genes ####
options(ggrepel.max.overlaps = 20)
cortex.integrated.sample.markers$log10padj <- -log10(cortex.integrated.sample.markers$p_val_adj)
cortex.integrated.sample.markers$avg_expr <- (cortex.integrated.sample.markers$pct.1 + cortex.integrated.sample.markers$pct.2)/2

for (i in c("MFG.03.41", "MFG.10.26", "MFG.10.62", "MFG.97.40")) {
  pdf(paste0("cortex.integrated.maplot.sample.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.sample.markers[cortex.integrated.sample.markers$cluster==i, ],
               mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=subset(cortex.integrated.sample.markers[cortex.integrated.sample.markers$cluster==i, ],
                                      (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                          aes(avg_expr, avg_log2FC, label=gene)) +
          ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Included") +
          xlab("Average Percent of Spots Expressing Gene") +
          ylab(paste0("Average Fold Change Between ", i, " and Others")))
  dev.off()
}

cortex.integrated.sample.markers.nomito <- cortex.integrated.sample.markers[!grepl('^MT-', cortex.integrated.sample.markers$gene),]

for (i in c("MFG.03.41", "MFG.10.26", "MFG.10.62", "MFG.97.40")) {
  pdf(paste0("cortex.integrated.maplot.sample.nomito.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.sample.markers.nomito[cortex.integrated.sample.markers.nomito$cluster==i, ],
               mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=subset(cortex.integrated.sample.markers.nomito[cortex.integrated.sample.markers.nomito$cluster==i, ],
                                      (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                          aes(avg_expr, avg_log2FC, label=gene)) +
          ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Excluded") +
          xlab("Average Percent of Spots Expressing Gene") +
          ylab(paste0("Average Fold Change Between ", i, " and Others")))
  dev.off()
}

#### MA plots by cluster instead of sample ####
cortex.integrated.cluster.markers$log10padj <- -log10(cortex.integrated.cluster.markers$p_val_adj)
cortex.integrated.cluster.markers$avg_expr <- (cortex.integrated.cluster.markers$pct.1 + cortex.integrated.cluster.markers$pct.2)/2
cortex.integrated.cluster.markers.nomito <- cortex.integrated.cluster.markers[!grepl('^MT-', cortex.integrated.cluster.markers$gene),]

for (i in c(0:13)) {
  pdf(paste0("cortex.integrated.maplot.cluster.nomito.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.cluster.markers.nomito[cortex.integrated.cluster.markers.nomito$cluster==i, ],
         mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
    geom_point() +
    scale_color_manual(name = 'p-adj < 1e-10',
                       values = setNames(c('blue', 'black'), c(T, F))) +
    geom_text_repel(data=subset(cortex.integrated.cluster.markers.nomito[cortex.integrated.cluster.markers.nomito$cluster==i, ],
                                (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                    aes(avg_expr, avg_log2FC, label=gene)) +
    ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Excluded") +
    xlab("Average Percent of Spots Expressing Gene") +
    ylab(paste0("Average Fold Change Between Cluster", i, " and Others")))
  dev.off()
}

for (i in c(0:13)) {
  pdf(paste0("cortex.integrated.maplot.cluster.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = cortex.integrated.cluster.markers[cortex.integrated.cluster.markers$cluster==i, ],
               mapping = aes(x = avg_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=subset(cortex.integrated.cluster.markers[cortex.integrated.cluster.markers$cluster==i, ],
                                      (avg_log2FC > 0.8 | avg_log2FC < -1) & avg_expr > 0.2),
                          aes(avg_expr, avg_log2FC, label=gene)) +
          ggtitle("Average Log2 Fold Change vs. Average Percent Spot Expression", subtitle="Mitochondrial Features Included") +
          xlab("Average Percent of Spots Expressing Gene") +
          ylab(paste0("Average Fold Change Between Cluster", i, " and Others")))
  dev.off()
}

# write feature expression tables
write.csv(cortex.integrated.cluster.markers, file="cortex.integrated.feature.expression.by.cluster.csv")
write.csv(cortex.integrated.sample.markers, file="cortex.integrated.feature.expression.by.sample.csv")
