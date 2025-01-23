#### LOAD LIBRARIES ####
pkg1 <- c('pheatmap', 'viridis', 'DESeq2', 'RColorBrewer', 'biomaRt', 'ggrepel','dplyr', 'optparse', 'ggplot2', 'tibble','stringr')
lapply(pkg1,
       function (x) library(x, character.only = TRUE))

#### READ IN ARGUMENTS ####
option_list <- list(
  make_option(c("-d", "--directory"), type="character",
              default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_6/short-DEGs-2",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character",
              default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_6/short-stringtie/gene_count_matrix.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_6/comparisons.csv",
              help="comparisons for differential expression")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$directory)

gcnts <-read.csv(opts$genes, header=TRUE, row.names = 1, check.names=FALSE)
print(colnames(gcnts))
#gcnts %>%
#  select(sort(names(.)))

gcnts <- gcnts[ , order(colnames(gcnts)) ]
print(colnames(gcnts))

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = 1)

comparisons <- comparisons[ order(row.names(comparisons)), ]
print(rownames(comparisons))

comparisons[comparisons == -1] <- 2

# standardize names (comparison names should be the full name, not the lab sample ID)
#colnames(gcnts) <- rownames(comparisons)

factors <- lapply(comparisons,
                  factor)

for (i in names(factors)) {
  comparisons[[i]] <- factor(comparisons[[i]])
}

dds <- sapply(names(comparisons),
              function (x) {
                f <- paste("~",x)
                print(f)
                DESeqDataSetFromMatrix(countData = gcnts,
                                       colData = comparisons,
                                       design = as.formula(f))
              },
              simplify = FALSE,
              USE.NAMES = TRUE
)

dds <- sapply(names(dds),
              function (x) {
                keep <- rowSums(counts(dds[[x]])) >= length(gcnts)
                return(dds[[x]][keep,])
              },
              simplify = FALSE,
              USE.NAMES = TRUE)

dds <- lapply(dds,
              DESeq)

save(dds, file = "deseq_dds_results.rData")


res <- sapply(names(dds),
              function (x) {
                results(dds[[x]], contrast = c(x,2,1))
              },
              simplify = FALSE,
              USE.NAMES = TRUE)

save(res, file = "deseq_res_results.rData")

resOrdered <- lapply(res,
                     function (x) {
                       x[order(x$padj),]
                       })

resSig <- lapply(resOrdered,
                 function (x) {
                   subset(x, padj < 0.05)
                 })

deNamesSig <- lapply(resSig,
                     rownames)

norm_counts <- sapply(names(dds),
                      function (x) {
                        counts(dds[[x]], normalized = TRUE)
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

write.table(norm_counts[[1]], file="deseq2_normalized_counts.txt", sep="\t", quote=F, col.names=NA)

normCountsSig <- sapply(names(norm_counts),
                        function (x) {
                          subset(norm_counts[[x]], rownames(norm_counts[[x]]) %in% deNamesSig[[x]])
                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE)

sapply(names(res),
       function (x) {
         ptNum=min(20,length(resSig[[x]]$baseMean))
         pdf(paste0(x, ".deseq.MAplot.pdf"))
         plotMA(res[[x]], xlim=c(1,1e5),
                ylim=c(-10,10), main=paste0(x), colSig="grey", colNonSig="grey")
         with(subset(res[[x]],
                     rownames(res[[x]]) %in% rownames(resSig[[x]])), {
          points(resSig[[x]]$baseMean, resSig[[x]]$log2FoldChange,
                 col=ifelse(resSig[[x]]$log2FoldChange>0,"forestgreen","red"),
                 cex = 0.4, lwd = 1, pch=19)
         })
         abline(h = c(-1,1), col = "blue")
         dev.off()
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
         results.ordered[[x]] %>%
           ggplot(aes(x = results.ordered[[x]]$log2FoldChange,
                      y = -log10(results.ordered[[x]]$padj))) +
           geom_point() + geom_abline(intercept = 2, slope = 0) +
           xlab("log2FoldChange") + ylab("-log10 padj") + ggtitle(paste0(x))
         ggsave(paste0(x, ".deseq.volcano.png"))
         })

#typically use rlog instead of vst, unless rlog throws a warning
ddst <- lapply(dds,
               function (x) {
                 rlog(x, blind = FALSE)
               })

save(ddst, file = "deseq_ddst_results.rData")

y <- names(comparisons)[1]
y <- str_split(y, '_v_')
contrast <- y[[1]][1]
base <- y[[1]][2]

sapply(names(ddst),
       function (x) {
         pdf(paste0(x,".DEseq2_PCA.pdf"))
         y <- x[1]
         y <- str_split(y, '_v_')
         y <- as.data.frame(y)
         contrast <- y[[1]][1]
         base <- y[[1]][2]
         pd <- DESeq2::plotPCA(ddst[[x]], intgroup = toString(x)) +
           scale_color_hue(labels = c("other", base, contrast))
         print(pd)
         dev.off()
       })

sampleDists <- lapply(ddst,
                      function (y) {
                        dist(t(assay(y)))
                      })
sampleDistMatrix <- lapply(sampleDists,
                           as.matrix)

topVarGenes <- head(order(rowVars(assay(dds[[1]])), decreasing = TRUE), 50)
minidat  <- assay(dds[[1]])[ topVarGenes, ]
pdf("topgenes.deseq.heatmap.pdf")
pheatmap(minidat,
         scale="row",
         color=inferno(30),
         fontsize_row = 4,
         fontsize_col = 6,
         main=paste0("Normalized Counts for Top 50 Variable Genes"),
)
dev.off()

sapply(names(resSig),
       function (x) {
         data <- as.data.frame(resSig[[x]])
         down <- subset(data, data[["log2FoldChange"]] < 0)
         up <- subset(data, data[["log2FoldChange"]] > 0)
         write.table(down, file=paste0(x,".deg.down.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(up, file=paste0(x,".deg.up.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
         write.table(data, file=paste0(x,".deg.all.sig.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })

sapply(names(res),
       function (x) {
         data <- as.data.frame(res[[x]])
         write.table(data, file=paste0(x, ".deg.allGenes.deseq.txt"), quote=FALSE, sep="\t", row.names = TRUE, col.names = NA)
       })
