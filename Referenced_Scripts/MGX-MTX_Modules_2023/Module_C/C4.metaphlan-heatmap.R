pkg1 <- c('pheatmap', 'viridis', 'RColorBrewer', 'biomaRt', 'dplyr', 'ggplot2', 'reshape2', 'cluster', 'gplots', 'tidyverse')
lapply(pkg1,
       function (x) library(x, character.only = TRUE))

newNames <- c("MB-001", "MB-003", "MB-004", "MB-005", "MB-006", "MB-011", "MB-012",
              "MB-014", "MB-016", "MB-017", "MB-018", "MB-020", "MB-021", "MB-023",
              "MB-024", "MB-025", "MB-028", "MB-033", "MB-035", "MB-037", "MB-038",
              "MB-039", "MB-040", "MB-044", "MB-045", "MB-047", "MB-049", "MB-050",
              "MB-053", "MB-055", "variance", "stdev")
rowMD <- c(rep("blue", 6), rep("red", 1),
           rep("blue", 4), rep("red", 1),
           rep("blue", 6), rep("red", 2),
           rep("blue", 10))

newNames <- c("49", "52", "56", "58", "60", "62", "variance", "stdev")
rowMD <- c(rep("blue", 2), rep("red", 2), rep("green", 2))

#### Set Variables ####

setwd("/Volumes/Gencore/analysis_projects/6078853_Otak_RNA/assembly-free-taxonomy/humann-prediction-pipeline/humann3_uniref50/taxonomic-classifications")

#### HEATMAP with full metaphlan dataset at SPECIES level ####

metaphlanS <- read.delim("merged_metaphlan_table_species.txt", header = TRUE, sep = "\t", row.names=1)
samplenum <- length(colnames(metaphlanS))
metaphlanS <- metaphlanS[rowSums(metaphlanS[])>0,]
metaphlanS$variance <- apply(metaphlanS,1,var)
metaphlanS$stdev <- apply(metaphlanS[,1:samplenum],1,sd, na.rm=TRUE)
metaphlanS <- metaphlanS[metaphlanS$stdev!=0,]
colnames(metaphlanS) <- newNames

hc <- hclust(as.dist(1-cor(metaphlanS[,1:samplenum], method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
pdf("sampleClusters-species.pdf")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")
dev.off()

hr <- hclust(as.dist(1-cor(t(metaphlanS[,1:samplenum]), method="pearson")), method="complete")
geneTree = as.dendrogram(hr, method="average")
pdf("speciesClusters.pdf")
plot(geneTree,
     leaflab = "none",             
     main = "Species Clustering",
     ylab = "Height")
dev.off()

metaphlanS<-as.matrix(metaphlanS)

wss <- (nrow(metaphlanS[,1:samplenum])-1)*sum(apply(metaphlanS[,1:samplenum],2,var))
for (i in 2:30) wss[i] <- sum(kmeans(metaphlanS[,1:samplenum],
                                     centers=i)$withinss)
pdf("sumSquares-species.pdf")
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

clusters <- 4

hclustk = cutree(hr, k=clusters)
clustColBar <- rainbow(length(unique(hclustk)), start=0.1, end=0.9)
clustColBar <- clustColBar[as.vector(hclustk)]

pdf(paste0("metaphlan-species.",clusters,".ClusteredHeatmap.pdf"), height=7, width=7)
heatmap.2(metaphlanS[,1:samplenum],
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=viridis(100),
          scale="row",
          cexCol = 1,
          labRow = F,
          main = "Metaphlan Clustered Species",
          trace = "none",
          ColSideColors=rowMD,
          RowSideColors=clustColBar,
          key = FALSE)
dev.off()

metaphlanS <- as.data.frame(metaphlanS)
metaphlanS$species <- rownames(metaphlanS)
metaphlanS$cluster <- hclustk
rownames(metaphlanS) <- NULL
row.order <- order.dendrogram(as.dendrogram(hr))
metaphlan.Ordered <- metaphlanS[match(row.order, rownames(metaphlanS)),]
metaphlan.Rev <- metaphlan.Ordered[nrow(metaphlan.Ordered):1,]
write.table(metaphlan.Rev, file="metaphlan_species_clusters.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

#### HEATMAP with full metaphlan dataset at GENUS level ####

metaphlanG <- read.delim("merged_metaphlan_table_genus.txt", header = TRUE, sep = "\t", row.names=1)
samplenum <- length(colnames(metaphlanG))
metaphlanG <- metaphlanG[rowSums(metaphlanG[])>0,]
metaphlanG$variance <- apply(metaphlanG,1,var)
metaphlanG$stdev <- apply(metaphlanG[,1:samplenum],1,sd, na.rm=TRUE)
metaphlanG <- metaphlanG[metaphlanG$stdev!=0,]
colnames(metaphlanG) <- newNames

hc <- hclust(as.dist(1-cor(metaphlanG[,1:samplenum], method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
pdf("sampleClusters-genus.pdf")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")
dev.off()

hr <- hclust(as.dist(1-cor(t(metaphlanG[,1:samplenum]), method="pearson")), method="complete")
geneTree = as.dendrogram(hr, method="average")
pdf("genusClusters.pdf")
plot(geneTree,
     leaflab = "none",             
     main = "genus Clustering",
     ylab = "Height")
dev.off()

metaphlanG<-as.matrix(metaphlanG)

wss <- (nrow(metaphlanG[,1:samplenum])-1)*sum(apply(metaphlanG[,1:samplenum],2,var))
for (i in 2:30) wss[i] <- sum(kmeans(metaphlanG[,1:samplenum],
                                     centers=i)$withinss)
pdf("sumSquares-genus.pdf")
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

clusters <- 24

hclustk = cutree(hr, k=clusters)
clustColBar <- rainbow(length(unique(hclustk)), start=0.1, end=0.9)
clustColBar <- clustColBar[as.vector(hclustk)]

pdf(paste0("metaphlan-genus.",clusters,".ClusteredHeatmap.pdf"), height=7, width=7)
heatmap.2(metaphlanG[,1:samplenum],
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=viridis(100),
          scale="row",
          cexCol = 1,
          labRow = F,
          main = "Metaphlan Clustered genus",
          trace = "none",
          ColSideColors=rowMD,
          RowSideColors=clustColBar,
          key = FALSE)
dev.off()

metaphlanG <- as.data.frame(metaphlanG)
metaphlanG$genus <- rownames(metaphlanG)
metaphlanG$cluster <- hclustk
rownames(metaphlanG) <- NULL
row.order <- order.dendrogram(as.dendrogram(hr))
metaphlan.Ordered <- metaphlanG[match(row.order, rownames(metaphlanG)),]
metaphlan.Rev <- metaphlan.Ordered[nrow(metaphlan.Ordered):1,]
write.table(metaphlan.Rev, file="metaphlan_genus_clusters.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

#### HEATMAP with full metaphlan dataset at FAMILY level ####

metaphlanF <- read.delim("merged_metaphlan_table_family.txt", header = TRUE, sep = "\t", row.names=1)
samplenum <- length(colnames(metaphlanF))
metaphlanF <- metaphlanF[rowSums(metaphlanF[])>0,]
metaphlanF$variance <- apply(metaphlanF,1,var)
metaphlanF$stdev <- apply(metaphlanF[,1:samplenum],1,sd, na.rm=TRUE)
metaphlanF <- metaphlanF[metaphlanF$stdev!=0,]
colnames(metaphlanF) <- newNames

hc <- hclust(as.dist(1-cor(metaphlanF[,1:samplenum], method="spearman")), method="complete")
sampleTree = as.dendrogram(hc, method="average")
pdf("sampleClusters-family.pdf")
plot(sampleTree,
     main = "Sample Clustering",
     ylab = "Height")
dev.off()

hr <- hclust(as.dist(1-cor(t(metaphlanF[,1:samplenum]), method="pearson")), method="complete")
geneTree = as.dendrogram(hr, method="average")
pdf("familyClusters.pdf")
plot(geneTree,
     leaflab = "none",             
     main = "family Clustering",
     ylab = "Height")
dev.off()

metaphlanF<-as.matrix(metaphlanF)

wss <- (nrow(metaphlanF[,1:samplenum])-1)*sum(apply(metaphlanF[,1:samplenum],2,var))
for (i in 2:30) wss[i] <- sum(kmeans(metaphlanF[,1:samplenum],
                                     centers=i)$withinss)
pdf("sumSquares-family.pdf")
plot(1:30, wss, type="b", xlab="Number of Clusters",
     ylab="Within groups sum of squares")
dev.off()

clusters <- 10

hclustk = cutree(hr, k=clusters)
clustColBar <- rainbow(length(unique(hclustk)), start=0.1, end=0.9)
clustColBar <- clustColBar[as.vector(hclustk)]

pdf(paste0("metaphlan-family.",clusters,".ClusteredHeatmap.pdf"), height=7, width=7)
heatmap.2(metaphlanF[,1:samplenum],
          Rowv=as.dendrogram(hr), 
          Colv=as.dendrogram(hc),
          col=viridis(100),
          scale="row",
          cexCol = 1,
          labRow = F,
          main = "Metaphlan Clustered family",
          trace = "none",
          ColSideColors=rowMD,
          RowSideColors=clustColBar,
          key = FALSE)
dev.off()

metaphlanF <- as.data.frame(metaphlanF)
metaphlanF$family <- rownames(metaphlanF)
metaphlanF$cluster <- hclustk
rownames(metaphlanF) <- NULL
row.order <- order.dendrogram(as.dendrogram(hr))
metaphlan.Ordered <- metaphlanF[match(row.order, rownames(metaphlanF)),]
metaphlan.Rev <- metaphlan.Ordered[nrow(metaphlan.Ordered):1,]
write.table(metaphlan.Rev, file="metaphlan_family_clusters.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
