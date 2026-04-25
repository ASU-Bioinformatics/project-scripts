library(pheatmap)
library(DESeq2)
library(dplyr)
library(viridis)
library(matrixStats)

setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp1/check")

# experimental DEseq2 conditions
gcnts <-read.csv("/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp1/gene_count_matrix.csv",
                 header=TRUE, row.names = 1, check.names=FALSE)

conditions <- read.csv("exp1-conditions.csv", header=TRUE, row.names = NULL)
print(conditions)

conditions$Treatment <- factor(conditions$Treatment, )
conditions$Day <- factor(conditions$Day)

conditions$Treatment <-relevel(conditions$Treatment, "Saline")

dds <- DESeqDataSetFromMatrix(countData = gcnts,
                                         colData = conditions,
                                         design = ~ Treatment + Day + Treatment:Day + rRNA)

dds <- estimateSizeFactors(dds)
nc <- counts(dds, normalized=TRUE)
filter <- rowSums(nc >= 10) >= 6
dds <- dds[filter,]

dds.deseq <- estimateSizeFactors(dds)
dds.deseq <- estimateDispersions(dds.deseq)
dds.deseq <- nbinomWaldTest(dds.deseq, maxit=1000)

ddsClean <- dds.deseq[which(mcols(dds.deseq)$betaConv),]

resultsNames(ddsClean)  

res.list <- c(
  "res.silk.v.saline.2" = results(ddsClean, name = "Treatment_Silk_vs_Saline"),
  "res.h1.v.saline.2" = results(ddsClean, name = "Treatment_H1.ag_vs_Saline"),
  "res.h4.v.saline.2" = results(ddsClean, name = "Treatment_H4.ag_vs_Saline"),
  "res.silkH1.v.saline.2" = results(ddsClean, name = "Treatment_Silk.plus.H1.ag_vs_Saline"),
  "res.silkH4.v.saline.2" = results(ddsClean, name = "Treatment_Silk.plus.H4.ag_vs_Saline"),
  "res.silk.v.saline.11" = results(ddsClean, contrast = list(c("Treatment_Silk_vs_Saline", "TreatmentSilk.Day11"))),
  "res.h1.v.saline.11" = results(ddsClean, contrast = list(c("Treatment_H1.ag_vs_Saline", "TreatmentH1.ag.Day11"))),
  "res.h4.v.saline.11" = results(ddsClean, contrast = list(c("Treatment_H4.ag_vs_Saline", "TreatmentH4.ag.Day11"))),
  "res.silkH1.v.saline.11" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H1.ag_vs_Saline", "TreatmentSilk.plus.H1.ag.Day11"))),
  "res.silkH4.v.saline.11" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H4.ag_vs_Saline", "TreatmentSilk.plus.H4.ag.Day11"))),
  "res.11.v.2.saline" = results(ddsClean, name = "Day_11_vs_2"),
  "res.11.v.2.silk" = results(ddsClean, contrast = list(c("Day_11_vs_2", "TreatmentSilk.Day11"))),
  "res.11.v.2.h1" = results(ddsClean, contrast = list(c("Day_11_vs_2", "TreatmentH1.ag.Day11"))),
  "res.11.v.2.h4" = results(ddsClean, contrast = list(c("Day_11_vs_2", "TreatmentH4.ag.Day11"))),
  "res.11.v.2.silkH1" = results(ddsClean, contrast = list(c("Day_11_vs_2", "TreatmentSilk.plus.H1.ag.Day11"))),
  "res.11.v.2.silkH4" = results(ddsClean, contrast = list(c("Day_11_vs_2", "TreatmentSilk.plus.H4.ag.Day11"))),
  "res.h1.v.silk.2" = results(ddsClean, contrast = list(c("Treatment_H1.ag_vs_Saline", "Treatment_Silk_vs_Saline"))),
  "res.h4.v.silk.2" = results(ddsClean, contrast = list(c("Treatment_H4.ag_vs_Saline", "Treatment_Silk_vs_Saline"))),
  "res.silkH1.v.silk.2" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H1.ag_vs_Saline", "Treatment_Silk_vs_Saline"))),
  "res.silkH4.v.silk.2" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H4.ag_vs_Saline", "Treatment_Silk_vs_Saline"))),
  "res.h1.v.silk.11" = results(ddsClean, contrast = list(c("Treatment_H1.ag_vs_Saline", "Treatment_Silk_vs_Saline",
                                                           "TreatmentH1.ag.Day11", "TreatmentSilk.Day11"))),
  "res.h4.v.silk.11" = results(ddsClean, contrast = list(c("Treatment_H4.ag_vs_Saline", "Treatment_Silk_vs_Saline",
                                                           "TreatmentH4.ag.Day11", "TreatmentSilk.Day11"))),
  "res.silkH1.v.silk.11" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H1.ag_vs_Saline", "Treatment_Silk_vs_Saline",
                                                               "TreatmentSilk.plus.H1.ag.Day11", "TreatmentSilk.Day11"))),
  "res.silkH4.v.silk.11" = results(ddsClean, contrast = list(c("Treatment_Silk.plus.H4.ag_vs_Saline", "Treatment_Silk_vs_Saline",
                                                               "TreatmentSilk.plus.H4.ag.Day11", "TreatmentSilk.Day11")))
  
)

file.names <- c(
  "res.silk.v.saline.2" = "Silk vs Saline, Day 2",
  "res.h1.v.saline.2" = "H1.ag vs Saline, Day 2",
  "res.h4.v.saline.2" = "H4.ag vs Saline, Day 2",
  "res.silkH1.v.saline.2" = "Silk + H1.ag vs Saline, Day 2",
  "res.silkH4.v.saline.2" = "Silk + H4.ag vs Saline, Day 2",
  "res.silk.v.saline.11" = "Silk vs Saline, Day 11",
  "res.h1.v.saline.11" = "H1.ag vs Saline, Day 11",
  "res.h4.v.saline.11" = "H4.ag vs Saline, Day 11",
  "res.silkH1.v.saline.11" = "Silk + H1.ag vs Saline, Day 11",
  "res.silkH4.v.saline.11" = "Silk + H4.ag vs Saline, Day 11",
  "res.11.v.2.saline" = "Day 11 vs Day 2, Saline",
  "res.11.v.2.silk" = "Day 11 vs Day 2, Silk",
  "res.11.v.2.h1" = "Day 11 vs Day 2, H1.ag",
  "res.11.v.2.h4" = "Day 11 vs Day 2, H4.ag",
  "res.11.v.2.silkH1" = "Day 11 vs Day 2, Silk + H1.ag",
  "res.11.v.2.silkH4" = "Day 11 vs Day 2, Silk + H4.ag",
  "res.h1.v.silk.2" = "H1.ag vs Silk, Day 2",
  "res.h4.v.silk.2" = "H4.ag vs Silk, Day 2",
  "res.silkH1.v.silk.2" = "Silk + H1.ag vs Silk, Day 2",
  "res.silkH4.v.silk.2" = "Silk + H4.ag vs Silk, Day 2",
  "res.h1.v.silk.11" = "H1.ag vs Silk, Day 11",
  "res.h4.v.silk.11" = "H4.ag vs Silk, Day 11",
  "res.silkH1.v.silk.11" = "Silk + H1.ag vs Silk, Day 11",
  "res.silkH4.v.silk.11" = "Silk + H4.ag vs Silk, Day 11"
)

resOrdered <- lapply(res.list,
                     function (x) {
                       x[order(x$padj),]
                     })

resSig <- lapply(resOrdered,
                 function (x) {
                   subset(x, padj < 0.05 & abs(log2FoldChange) > 2)
                 })

deNamesSig <- lapply(resSig,
                     rownames)

norm_counts <- counts(ddsClean, normalized = TRUE)

write.table(norm_counts, file="deseq2_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

normCountsSig <- sapply(names(deNamesSig),
                        function (x) {
                          subset(norm_counts, rownames(norm_counts) %in% deNamesSig[[x]])
                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE)

sapply(names(res.list),
       function (x) {
         ptNum=min(20,length(resSig[[x]]$baseMean))
         pdf(paste0(x, ".MAplot.pdf"), height = 5, width = 7)
         DESeq2::plotMA(res.list[[x]], xlim=c(1,1e5),
                        ylim=c(-10,10), main=paste0(x), colSig="grey", colNonSig="grey")
         with(subset(res[[x]],
                     rownames(res.list[[x]]) %in% rownames(resSig[[x]])), {
                       points(resSig[[x]]$baseMean, resSig[[x]]$log2FoldChange,
                              col=ifelse(resSig[[x]]$log2FoldChange>2,"forestgreen","red"),
                              cex = 0.4, lwd = 1, pch=19)
                     })
         abline(h = c(-2,2), col = "blue")
         dev.off()
         ggplot(res.list[[x]], aes(baseMean, log2FoldChange, colour=padj)) + 
           geom_point(size=1) + 
           scale_y_continuous(limits=c(-6, 6), oob=squish) + 
           scale_x_log10() + 
           geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
           labs(x="mean of normalized counts", y="log fold change") + 
           scale_colour_viridis(direction=-1, trans = "sqrt") + 
           theme_bw() + 
           geom_density_2d(colour="black", size=1.5, bins=6) +
           ggtitle(file.names[[x]])
         ggsave(paste0(x, ".contouredMAplot.pdf"), height = 5, width = 7)
       })

results.ordered <- sapply(names(resOrdered),
                          function (x) {
                            x <- as.data.frame(resOrdered[[x]]) %>%
                              arrange(padj)
                          },
                          simplify = FALSE,
                          USE.NAMES = TRUE)


lapply(names(results.ordered),
       function (x) {
         results.ordered[[x]]$names <- row.names(results.ordered[[x]])
         results.ordered[[x]]$names <- results.ordered[[x]]$names <- (separate_wider_delim(results.ordered[[x]], 
                                                                                           "names", "|", 
                                                                                           names = c("id", "name"), 
                                                                                           too_few = "align_start"))$name
         results.ordered[[x]]$names[16:nrow(results.ordered[[x]])] <- ""
         results.ordered[[x]]$colors <- "grey"
         results.ordered[[x]]$colors[results.ordered[[x]]$padj<=0.05 & results.ordered[[x]]$log2FoldChange>2] <- "forestgreen"
         results.ordered[[x]]$colors[results.ordered[[x]]$padj<=0.05 & results.ordered[[x]]$log2FoldChange<2] <- "red"
         results.ordered[[x]] %>%
           ggplot(aes(x = results.ordered[[x]]$log2FoldChange,
                      y = -log10(results.ordered[[x]]$padj))) +
           geom_point(aes(colour = colors )) + scale_color_identity() + geom_abline(intercept = 2, slope = 0) +
           ggrepel::geom_text_repel(aes(label = names), 
                                    size = 3.0,
                                    point.padding = 0.1,
                                    box.padding = 0.75) + 
           theme_light() +
           xlab("log2FoldChange") + ylab("-log10 padj") + ggtitle(paste0(file.names[[x]]))
         ggsave(paste0(x, ".volcano.png"), width = 7, height = 7, units = "in", dpi = 720)
       })

ddsTransform <- vst(ddsClean, blind = FALSE)

sizeCol = ddsTransform@colData@listData[["rRNA"]]
shapeCol = ddsTransform@colData@listData[["Day"]]

DESeq2::plotPCA(ddsTransform, intgroup = "Treatment") +
  #geom_point(size=0.5) +
  geom_point(aes(shape = shapeCol, size=sizeCol)) +
  scale_color_discrete(palette = "Dark2") +
  scale_radius(range=c(5,25)) +
  theme_minimal(base_size = 18) +
  labs( title = "Principal Components",
        subtitle = "based on top 500 most variable genes, by transformed normalized counts",
        size = "rRNA abundance",
        shape = "day",
        color = "treatment")

ggsave("principal.components.png", width = 16, height = 12, units = "in", dpi = 720)

ddsVST <- as.data.frame(assay(ddsTransform))
ddsVST$Gene <- rownames(ddsVST)

sigGenes <- unique(c(rownames(subset(res.list$res.silk.v.saline.2, 
                                   padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h1.v.saline.2, 
                                   padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h4.v.saline.2, 
                                   padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH1.v.saline.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH4.v.saline.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silk.v.saline.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h1.v.saline.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h4.v.saline.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH1.v.saline.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH4.v.saline.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.saline, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.silk, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.h1, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.h4, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.silkH1, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.11.v.2.silkH4, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h1.v.silk.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h4.v.silk.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH1.v.silk.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH4.v.silk.2, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h1.v.silk.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.h4.v.silk.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH1.v.silk.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2)),
                     rownames(subset(res.list$res.silkH4.v.silk.11, 
                                     padj <= 0.05 & abs(log2FoldChange) >= 2))
))

minidat  <- assay(ddsTransform)[ sigGenes, ]

pdf("siggenes.heatmap.pdf")
pheatmap(minidat,
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         show_rownames = F,
         main=paste0("Normalized Counts for Significantly Differential Genes Across All Comparisons"),
)
dev.off()

minidat2 <- as.data.frame(minidat)
minidat2$var <- apply(minidat2, 1, var)
minidat2 <- minidat2[order(minidat2$var), ]

pdf("variable-siggenes.heatmap.pdf")
pheatmap(minidat2[c(1:50), c(1:48)],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         fontsize_row = 4,
         show_rownames = T,
         main=paste0("Normalized Counts for 50 Most Variable Significantly Differential Genes"),
)
dev.off()

sapply(names(resSig),
       function (x) {
         data <- as.data.frame(resSig[[x]])
         down <- subset(data, data[["log2FoldChange"]] < 0)
         up <- subset(data, data[["log2FoldChange"]] > 0)
         write.table(down, file=paste0(x,".deg.down.sig.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(up, file=paste0(x,".deg.up.sig.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(data, file=paste0(x,".deg.all.sig.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })

sapply(names(resOrdered),
       function (x) {
         data <- as.data.frame(resOrdered[[x]])
         write.table(data, file=paste0(x, ".deg.allGenes.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })
