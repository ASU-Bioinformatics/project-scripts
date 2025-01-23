# load module hdf5-1.13.1-gcc-11.2.0

install.packages("Seurat")
install.packages("readbitmap")
install.packages("hdf5r")
devtools::install_github('satijalab/seurat-data')
BiocManager::install('glmGamPoi')
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
MFG.03.41 <- Load10X_Spatial(data.dir = "/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Nafisa_MFG-03-41/outs")
MFG.10.26 <- Load10X_Spatial(data.dir = "/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Nafisa_MFG-10-26/outs")
MFG.10.62 <- Load10X_Spatial(data.dir = "/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Nafisa_MFG-10-62/outs")
MFG.97.40 <- Load10X_Spatial(data.dir = "/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Nafisa_MFG-97-40/outs")

#### MFG.10.26 ####
# the poor clustering here may be due to the larger percentage of genomic DNA (44.7%)

#### MFG.03.41 ####
plot1 <- VlnPlot(MFG.03.41, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(MFG.03.41, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

pdf("MFG-03-41.raw-counts.pdf")
wrap_plots(plot1, plot2)
dev.off()

MFG.03.41 <- PercentageFeatureSet(MFG.03.41, "^MT-", col.name = "percent_mito")
MFG.03.41 <- PercentageFeatureSet(MFG.03.41, "^RP", col.name = "percent_ribosomal")

pdf("MFG-03-41.general-stats-violin.pdf", width=12, height=12)
VlnPlot(MFG.03.41, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf("MFG-03-41.general-stats-spatial.pdf", width=16, height=12)
SpatialFeaturePlot(MFG.03.41, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"), ncol = 2)
dev.off()

MFG.03.41.filtered <- MFG.03.41[, MFG.03.41$nFeature_Spatial > 500 & MFG.03.41$percent_mito < 40]
pdf("MFG-03-41.filtered.general-stats-violin.pdf", width=12, height=12)
VlnPlot(MFG.03.41.filtered, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal"),
        pt.size = 0.1, ncol = 2) + NoLegend()
dev.off()

pdf("MFG-03-41.filtered.general-stats-spatial.pdf", width=16, height=12)
SpatialFeaturePlot(MFG.03.41.filtered, features = c("nCount_Spatial", "nFeature_Spatial", "percent_mito", "percent_ribosomal", ncol = 2))
dev.off()

C = MFG.03.41.filtered@assays$Spatial@layers$counts
rownames(C) <- Features(MFG.03.41.filtered)
C@x = C@x/rep.int(colSums(C), diff(C@p))
most_expressed <- order(Matrix::rowSums(C), decreasing = T)[20:1]

pdf("MFG-03-41.filtered.high-expression.pdf", width = 8, height = 5)
boxplot(as.matrix(t(C[most_expressed, ])), cex.axis = 0.6, las = 1, xlab = "% total count per cell",
        col = (scales::hue_pal())(20)[20:1], horizontal = TRUE, main="Most Expressed Genes Across MFG.03.041")
dev.off()

MFG.03.41.filtered <- SCTransform(MFG.03.41.filtered, assay = "Spatial", verbose = FALSE)

plot1 <- VlnPlot(MFG.03.41.filtered, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(MFG.03.41.filtered, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)

MFG.03.41.filtered <- RunPCA(MFG.03.41.filtered, assay = "SCT", verbose = FALSE)
MFG.03.41.filtered <- FindNeighbors(MFG.03.41.filtered, reduction = "pca", dims = 1:30)
MFG.03.41.filtered <- FindClusters(MFG.03.41.filtered, verbose = FALSE)
MFG.03.41.filtered <- RunUMAP(MFG.03.41.filtered, reduction = "pca", dims = 1:30)
MFG.03.41.filtered <- RunTSNE(MFG.03.41.filtered, reduction = "pca", dims = 1:30)

p1 <- DimPlot(MFG.03.41.filtered, reduction = "umap", label = TRUE, pt.size=0.1) + 
  NoLegend() + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9))
p2 <- DimPlot(MFG.03.41.filtered, reduction = "tsne", label = TRUE, pt.size=0.1) +
  NoLegend() + 
  theme(axis.text = element_text(size=8),
        axis.title = element_text(size=9))
p3 <- SpatialDimPlot(MFG.03.41.filtered, label = TRUE, label.size = 3)

pdf("MFG-03-41.filtered.normalized.clusters.pdf", width=12, height=4)
p3 + p1 + p2 + 
  plot_annotation(title ="Clusters Visualized as Spatial Overlay, UMAP Reduction, and tSNE Reduction") &
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

MFG.03.41.filtered <- FindSpatiallyVariableFeatures(MFG.03.41.filtered, assay = "SCT",
                                                    features = VariableFeatures(MFG.03.41.filtered)[1:1000],
                                                    selection.method = "moransi")
top.features <- head(SpatiallyVariableFeatures(MFG.03.41.filtered, selection.method = "moransi"), 12)

MFG.03.41.filtered <- FindSpatiallyVariableFeatures(MFG.03.41.filtered, assay = 'SCT',
                                                    features = VariableFeatures(MFG.03.41.filtered)[1:1000],
                                                    selection.method = 'markvariogram')
top.features.mv <- head(SpatiallyVariableFeatures(MFG.03.41.filtered, selection.method = "markvariogram"), 12)

SpatialColors <- colorRampPalette(colors = rev(x = brewer.pal(n = 11, name = "Spectral")))

pdf("MFG-03-41.filtered.top12features.moransi.pdf", width=20, height=20)
SpatialFeaturePlot(MFG.03.41.filtered, features = top.features, 
                   ncol = 4, alpha = c(0.5, 1), pt.size.factor = 3) & 
  ggplot2::scale_fill_gradientn(limits = c(0.0, 5.0), colours = SpatialColors(n=100))
dev.off()

pdf("MFG-03-41.filtered.top12features.markvariogram.pdf", width=20, height=20)
SpatialFeaturePlot(MFG.03.41.filtered, features = top.features.mv, 
                   ncol = 4, alpha = c(0.5, 1), pt.size.factor = 3) & 
  ggplot2::scale_fill_gradientn(limits = c(0.0, 5.0), colours = SpatialColors(n=100))
dev.off()

#### MFG.97.40 ####
