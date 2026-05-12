#### LOAD LIBRARIES ####
pkg1 <- c('pheatmap', 'viridis', 'DESeq2', 'RColorBrewer', 'biomaRt', 'ggrepel','dplyr', 'optparse', 'ggplot2', 'tibble','stringr', 'scales','tidyr')
lapply(pkg1,
       function (x) library(x, character.only = TRUE))

#### READ IN ARGUMENTS ####
option_list <- list(
  make_option(c("-d", "--directory"), type="character",
              default="/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp2",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character",
              default="/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp2/gene_count_matrix.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp2/exp2-comparisons.csv",
              help="comparisons for differential expression")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$directory)

gcnts <-read.csv(opts$genes, header=TRUE, row.names = 1, check.names=FALSE)
print(colnames(gcnts))

gcnts <- gcnts[ , order(colnames(gcnts)) ]
print(colnames(gcnts))

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = NULL)
print(comparisons)

comparisons <- comparisons[ order(comparisons$Name), ]
print(comparisons$Name)

compRowNames <- comparisons$Name
compColNames <- colnames(comparisons)
comparisons <- as.data.frame(comparisons[ , -1])
row.names(comparisons) <- compRowNames
colnames(comparisons) <- compColNames[-1]

comparisons[comparisons == -1] <- 2

factors <- lapply(comparisons,
                  factor)

for (i in names(factors)) {
  comparisons[[i]] <- factor(comparisons[[i]])
}

comparisons$rRNA <- c(0.009129917, 0.03727596, 0.008394049, 0.711475743, 0.661464187, 0.256237324) # this is total number of alignments from the GTF w/o rRNA divided by the total number of alignments from the GTF w/ rRNA
compnum <- length(comparisons) - 2

#### MODELING AND DEG rRNA MAIN EFFECT ####
dds.rRNAonly <- sapply(names(comparisons[1:compnum]),
                       function (x) {
                         f <- paste("~",x,"+rRNA")
                         print(f)
                         DESeqDataSetFromMatrix(countData = gcnts,
                                                colData = comparisons,
                                                design = as.formula(f))
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE
)

dds.rRNAonly <- sapply(names(dds.rRNAonly),
                       function (x) {
                         keep <- rowSums(counts(dds.rRNAonly[[x]]) > 5 ) >= length(gcnts)
                         return(dds.rRNAonly[[x]][keep,])
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE)

dds.rRNAonly <- lapply(dds.rRNAonly,
                       DESeq)

res.rRNAonly <- sapply(names(dds.rRNAonly),
                       function (x) {
                         results(dds.rRNAonly[[x]], contrast = c(x,2,1))
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE)

resOrdered.rRNAonly <- lapply(res.rRNAonly,
                              function (x) {
                                x[order(x$padj),]
                              })

resSig.rRNAonly <- lapply(resOrdered.rRNAonly,
                          function (x) {
                            subset(x, padj < 0.05)
                          })

deNamesSig.rRNAonly <- lapply(resSig.rRNAonly,
                              rownames)

norm_counts.rRNAonly <- sapply(names(dds.rRNAonly),
                               function (x) {
                                 counts(dds.rRNAonly[[x]], normalized = TRUE)
                               },
                               simplify = FALSE,
                               USE.NAMES = TRUE)

write.table(norm_counts.rRNAonly[[1]], file="deseq2_normalized_counts.rRNAonly.txt", sep="\t", quote=F, col.names=NA)

normCountsSig.rRNAonly <- sapply(names(norm_counts.rRNAonly),
                                 function (x) {
                                   subset(norm_counts.rRNAonly[[x]], rownames(norm_counts.rRNAonly[[x]]) %in% deNamesSig.rRNAonly[[x]])
                                 },
                                 simplify = FALSE,
                                 USE.NAMES = TRUE)

sapply(names(res.rRNAonly),
       function (x) {
         ptNum=min(20,length(resSig.rRNAonly[[x]]$baseMean))
         pdf(paste0(x, ".rRNAonly.MAplot.pdf"), height = 5, width = 7)
         DESeq2::plotMA(res.rRNAonly[[x]], xlim=c(1,1e5),
                        ylim=c(-10,10), main=paste0(x), colSig="grey", colNonSig="grey")
         with(subset(res.rRNAonly[[x]],
                     rownames(res.rRNAonly[[x]]) %in% rownames(resSig.rRNAonly[[x]])), {
                       points(resSig.rRNAonly[[x]]$baseMean, resSig.rRNAonly[[x]]$log2FoldChange,
                              col=ifelse(resSig.rRNAonly[[x]]$log2FoldChange>2,"forestgreen","red"),
                              cex = 0.4, lwd = 1, pch=19)
                     })
         abline(h = c(-2,2), col = "blue")
         dev.off()
         ggplot(res.rRNAonly[[x]], aes(baseMean, log2FoldChange, colour=padj)) + 
           geom_point(size=1) + 
           scale_y_continuous(limits=c(-6, 6), oob=squish) + 
           scale_x_log10() + 
           geom_hline(yintercept = 0, colour="darkorchid4", size=1, linetype="longdash") + 
           labs(x="mean of normalized counts", y="log fold change") + 
           scale_colour_viridis(direction=-1, trans = "sqrt") + 
           theme_bw() + 
           geom_density_2d(colour="black", size=1.5, bins=6) +
           ggtitle("MpDNA vs UMpDNA")
         ggsave(paste0(x, ".contouredMAplot.pdf"), height = 5, width = 7)
       })

results.ordered.rRNAonly <- sapply(names(resOrdered.rRNAonly),
                                   function (x) {
                                     x <- as.data.frame(resOrdered.rRNAonly[[x]]) %>%
                                       arrange(padj)
                                   },
                                   simplify = FALSE,
                                   USE.NAMES = TRUE)

lapply(names(results.ordered.rRNAonly),
       function (x) {
         results.ordered.rRNAonly[[x]]$names <- row.names(results.ordered.rRNAonly[[x]])
         results.ordered.rRNAonly[[x]]$names <- results.ordered.rRNAonly[[x]]$names <- (separate_wider_delim(results.ordered.rRNAonly[[x]], 
                                                                                                             "names", "|", 
                                                                                                             names = c("id", "name"), 
                                                                                                             too_few = "align_start"))$name
         results.ordered.rRNAonly[[x]]$names[16:nrow(results.ordered.rRNAonly[[x]])] <- ""
         results.ordered.rRNAonly[[x]]$colors <- "grey"
         results.ordered.rRNAonly[[x]]$colors[results.ordered.rRNAonly[[x]]$padj<=0.05 & results.ordered.rRNAonly[[x]]$log2FoldChange>2] <- "forestgreen"
         results.ordered.rRNAonly[[x]]$colors[results.ordered.rRNAonly[[x]]$padj<=0.05 & results.ordered.rRNAonly[[x]]$log2FoldChange<2] <- "red"
         results.ordered.rRNAonly[[x]] %>%
           ggplot(aes(x = results.ordered.rRNAonly[[x]]$log2FoldChange,
                      y = -log10(results.ordered.rRNAonly[[x]]$padj))) +
           geom_point(aes(colour = colors )) + scale_color_identity() + geom_abline(intercept = 2, slope = 0) +
           ggrepel::geom_text_repel(aes(label = names), 
                                    size = 3.0,
                                    point.padding = 0.1,
                                    box.padding = 0.75) + 
           theme_light() +
           xlab("log2FoldChange") + ylab("-log10 padj") + ggtitle(paste0("MpDNA vs UMpDNA"))
         ggsave(paste0(x, ".volcano.png"), width = 7, height = 7, units = "in", dpi = 720)
       })

ddsTransform <- vst(dds.rRNAonly$MpDNA_v_UMpDNA, blind = FALSE)

sizeCol = ddsTransform@colData@listData[["rRNA"]]

DESeq2::plotPCA(ddsTransform, intgroup = "Type") +
  #geom_point(size=0.5) +
  geom_point(aes(size=sizeCol)) +
  scale_color_discrete(palette = "Dark2") +
  scale_radius(range=c(5,25)) +
  theme_minimal(base_size = 18) +
  labs( title = "Principal Components",
        subtitle = "based on top 500 most variable genes, by transformed normalized counts",
        size = "rRNA abundance",
        color = "treatment")

ggsave("principal.components.png", width = 16, height = 12, units = "in", dpi = 720)

#typically use rlog instead of vst, unless rlog throws a warning
tryCatch( { dds.rRNAonlyt <- lapply(dds.rRNAonly,
                                    function (x) {
                                      rlog(x, blind = FALSE)
                                    }) } ,
          warning = function(w) { print("used vst")
            dds.rRNAonlyt <- lapply(dds.rRNAonly,
                                    function (x) {
                                      vst(x, blind = FALSE)
                                    })} )

save(dds.rRNAonlyt, file = "deseq_dds.rRNAonlyt_results.rData")

y <- names(comparisons)[1]
y <- str_split(y, '_v_')
contrast <- y[[1]][1]
base <- y[[1]][2]

sampleDists.rRNAonly <- lapply(dds.rRNAonlyt,
                               function (y) {
                                 dist(t(assay(y)))
                               })
sampleDistMatrix.rRNAonly <- lapply(sampleDists.rRNAonly,
                                    as.matrix)

topVarGenes.rRNAonly <- head(order(rowVars(assay(dds.rRNAonly[[1]])), decreasing = TRUE), 50)
minidat.rRNAonly  <- assay(dds.rRNAonly[[1]])[ topVarGenes.rRNAonly, ]
pdf("topgenes.rRNAonly.deseq.heatmap.pdf")
pheatmap(minidat.rRNAonly,
         scale="row",
         color=inferno(30),
         fontsize_row = 4,
         fontsize_col = 6,
         main=paste0("Normalized Counts for Top 50 Most Variable Genes"),
)
dev.off()

sapply(names(resSig.rRNAonly),
       function (x) {
         data <- as.data.frame(resSig.rRNAonly[[x]])
         down <- subset(data, data[["log2FoldChange"]] < 0)
         up <- subset(data, data[["log2FoldChange"]] > 0)
         write.table(down, file=paste0(x,".rRNAonly.deg.down.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(up, file=paste0(x,".rRNAonly.deg.up.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(data, file=paste0(x,".rRNAonly.deg.all.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })

sapply(names(res.rRNAonly),
       function (x) {
         data <- as.data.frame(res.rRNAonly[[x]])
         write.table(data, file=paste0(x, ".rRNAonly.deg.allGenes.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })
