# load module hdf5-1.13.1-gcc-11.2.0
library(Seurat)
library(SeuratData)
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

#### read in all samples and create Seurat Objects ####
MFG.97.40 <- Load10X_Spatial(data.dir = "/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Nafisa_MFG-97-40/outs")

#### MFG.97.40 ####
plot1 <- VlnPlot(MFG.97.40, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(MFG.97.40, features = "nCount_Spatial", pt.size.factor = 3) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

pdf("MFG-97-40.raw-counts.pdf")
wrap_plots(plot1, plot2)
dev.off()

MFG.97.40 <- PercentageFeatureSet(MFG.97.40, "^MT-", col.name = "percent_mito")
MFG.97.40 <- PercentageFeatureSet(MFG.97.40, "^RP", col.name = "percent_ribosomal")

pdf("MFG-97-40.general-stats-violin.pdf", width=12, height=12)
VlnPlot(MFG.97.40, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf("MFG-97-40.general-stats-spatial.pdf", width=16, height=12)
SpatialFeaturePlot(MFG.97.40,
                   features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
                   ncol = 2, pt.size.factor = 3)
dev.off()

# to recreate this plot after downstream analysis has been run, reset default assay and active idents
# Idents(MFG.97.40.filtered) <- MFG.97.40.filtered$orig.ident
# DefaultAssay(MFG.97.40.filtered) <- "Spatial"

# then, to reset them to the downstream state again
Idents(MFG.97.40.filtered) <- MFG.97.40.filtered$seurat_clusters
DefaultAssay(MFG.97.40.filtered) <- "SCT"

MFG.97.40.filtered <- MFG.97.40[, MFG.97.40$nFeature_Spatial > 500 & MFG.97.40$percent_mito < 40]
pdf("MFG-97-40.filtered.general-stats-violin.pdf", width=12, height=12)
VlnPlot(MFG.97.40.filtered, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend()  +
  plot_annotation(title="General Metrics for MFG-97-40", subtitle="Filtered to >500 features, <40% mitochondrial")
dev.off()

pdf("MFG-97-40.filtered.general-stats-spatial.pdf", width=16, height=12)
SpatialFeaturePlot(MFG.97.40.filtered,
                   features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
                   ncol = 2, pt.size.factor = 3) +
  plot_annotation(title="General Metrics for MFG-97-40", subtitle="Filtered to >500 features, <40% mitochondrial")
dev.off()

C = MFG.97.40.filtered@assays$Spatial@layers$counts
rownames(C) <- Features(MFG.97.40.filtered)
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

pdf("MFG-97-40.filtered.high-expression.pdf", width = 8, height = 5)
boxplot(as.matrix(t(C[most_expressed, ])), cex.axis = 0.6, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main="Most Expressed Genes Across MFG.97.40")
dev.off()

MFG.97.40.filtered <- SCTransform(MFG.97.40.filtered, assay = "Spatial", verbose = FALSE)

plot1 <- VlnPlot(MFG.97.40.filtered, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(MFG.97.40.filtered,
                            features = "nCount_Spatial",
                            pt.size.factor = 3) + theme(legend.position = "right")
wrap_plots(plot1, plot2)

MFG.97.40.filtered <- RunPCA(MFG.97.40.filtered, assay = "SCT", verbose = FALSE)
MFG.97.40.filtered <- FindNeighbors(MFG.97.40.filtered, reduction = "pca", dims = 1:30)
MFG.97.40.filtered <- FindClusters(MFG.97.40.filtered, verbose = FALSE)
MFG.97.40.filtered <- RunUMAP(MFG.97.40.filtered, reduction = "pca", dims = 1:30)
MFG.97.40.filtered <- RunTSNE(MFG.97.40.filtered, reduction = "pca", dims = 1:30)

p1 <- DimPlot(MFG.97.40.filtered, reduction = "umap", label = TRUE, pt.size=0.1) + 
  NoLegend() + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9))
p2 <- DimPlot(MFG.97.40.filtered, reduction = "tsne", label = TRUE, pt.size=0.1) +
  NoLegend() + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9))
p3 <- SpatialDimPlot(MFG.97.40.filtered, label = TRUE, label.size = 3, pt.size.factor = 3)

pdf("MFG-97-40.filtered.normalized.clusters.pdf", width=12, height=4)
p3 + p1 + p2 + 
  plot_annotation(title ="Clusters Visualized as Spatial Overlay, UMAP Reduction, and tSNE Reduction") &
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

MFG.97.40.filtered <- FindSpatiallyVariableFeatures(MFG.97.40.filtered, assay = "SCT",
                                                    features = VariableFeatures(MFG.97.40.filtered)[1:1000],
                                                    selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(MFG.97.40.filtered, selection.method = "moransi"), 12)

MFG.97.40.filtered <- FindSpatiallyVariableFeatures(MFG.97.40.filtered, assay = 'SCT',
                                                    features = VariableFeatures(MFG.97.40.filtered)[1:1000],
                                                    selection.method = 'markvariogram')
top.features.mv <- head(SpatiallyVariableFeatures(MFG.97.40.filtered, selection.method = "markvariogram"), 12)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

pdf("MFG-97-40.filtered.top12features.moransi.pdf", width=20, height=20)
SpatialFeaturePlot(MFG.97.40.filtered, features = top.features, 
                   ncol = 4, alpha = c(0.5, 1), pt.size.factor = 3)  +
  plot_annotation(title="Distribution for Top 12 Spatially Variable Features", subtitle="Moran's I Method, Excluding Mitochondrial Genes") & 
  ggplot2::scale_fill_gradientn(limits = c(0.0, 5.0), colours = SpatialColors(n=100))
dev.off()

pdf("MFG-97-40.filtered.top12features.markvariogram.pdf", width=20, height=20)
SpatialFeaturePlot(MFG.97.40.filtered, features = top.features.mv, 
                   ncol = 4, alpha = c(0.5, 1), pt.size.factor = 3) +
  plot_annotation(title="Distribution for Top 12 Spatially Variable Features", subtitle="Mark's Variogram Method, Excluding Mitochondrial Genes") & 
  ggplot2::scale_fill_gradientn(limits = c(0.0, 5.0), colours = SpatialColors(n=100))
dev.off()

#### find differential markers between clusters ####

MFG.97.40.filtered.markers <- FindAllMarkers(MFG.97.40.filtered, assay="SCT", )

MFG.97.40.filtered.markers$diff_expr <- abs(MFG.97.40.filtered.markers$pct.1 - MFG.97.40.filtered.markers$pct.2)
MFG.97.40.filtered.markers$mass <- abs(MFG.97.40.filtered.markers$diff_expr * MFG.97.40.filtered.markers$avg_log2FC)
MFG.97.40.filtered.markers.nomito <- MFG.97.40.filtered.markers[!grepl('^MT-', MFG.97.40.filtered.markers$gene),]

for (i in c(0:11)) {
  pdf(paste0("MFG.97.40.filtered.maplot.cluster.nomito.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = MFG.97.40.filtered.markers.nomito[MFG.97.40.filtered.markers.nomito$cluster==i, ],
               mapping = aes(x = diff_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=MFG.97.40.filtered.markers %>% 
                            arrange(-mass) %>%
                            filter(cluster==i) %>%
                            top_n(5),
                          aes(diff_expr, avg_log2FC, label=gene)) +
          ggtitle("MFG-97-40 Average Log2 Fold Change vs.Percent Spot Expression Difference", subtitle="Mitochondrial Features Excluded") +
          xlab("Difference between Percent of Spots Expressing Feature in MFG.97.40 and Others") +
          ylab(paste0("Average Fold Change Between Cluster", i, " and Others")))
  dev.off()
}

for (i in c(0:11)) {
  pdf(paste0("MFG.97.40.filtered.maplot.cluster.",i,".pdf"), width = 8, height = 6)
  print(ggplot(data = MFG.97.40.filtered.markers[MFG.97.40.filtered.markers$cluster==i, ],
               mapping = aes(x = diff_expr, y = avg_log2FC, color = p_val_adj < 0.0000000001)) +
          geom_point() +
          scale_color_manual(name = 'p-adj < 1e-10',
                             values = setNames(c('blue', 'black'), c(T, F))) +
          geom_text_repel(data=MFG.97.40.filtered.markers %>% 
                            arrange(-mass) %>%
                            filter(cluster==i) %>%
                            top_n(5),
                          aes(diff_expr, avg_log2FC, label=gene)) +
          ggtitle(" MFG-97-40 Average Log2 Fold Change vs. Percent Spot Expression Difference", subtitle="Mitochondrial Features Included") +
          xlab("Difference between Percent of Spots Expressing Feature in MFG.97.40 and Others") +
          ylab(paste0("Average Fold Change Between Cluster", i, " and Others")))
  dev.off()
}

# write feature expression tables
write.csv(MFG.97.40.filtered.markers, file="MFG.97.40.feature.expression.by.cluster.csv")
