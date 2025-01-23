#### LOAD LIBRARIES ####
pkg1 <- c('edgeR', 'viridis', 'pheatmap', 'optparse', 'dplyr','tibble')
lapply(pkg1,
       function (x) library(x, character.only = TRUE))

#### READ IN ARGUMENTS ####
option_list <- list(
  make_option(c("-d", "--directory"), type="character",
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/DEGS-lowqual",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character",
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual/gene_count_matrix.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/comparisons.csv",
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

comparisons <- comparisons[ order(comparisons$Name), ]
print(comparisons$Name)

compRowNames <- comparisons$Name
compColNames <- colnames(comparisons)
comparisons <- as.data.frame(comparisons[ , -1])
row.names(comparisons) <- compRowNames
colnames(comparisons) <- compColNames[-1]

#define groupings for edgeR
grouplist <- lapply(comparisons,
                    factor)

#define comparisons for edgeR
compList <- sapply(grouplist,
                   function (x) {
                     print(nlevels(x))
                     if (nlevels(x) == 3) {
                       x <- list(x, c(1,0,-1))
                     } else if (nlevels(x) == 2) {
                       x <- list(x, c(1,-1))
                     }
                   },
                   simplify = FALSE,
                   USE.NAMES = TRUE)

compList
head(gcnts)

#### EDGER: MODEL AND MDS PLOT ####
y<-DGEList(counts=gcnts)
y<-calcNormFactors(y)

num <- length(gcnts)

pal <- viridis_pal(option = "H", begin=0, end=0.9)(num)

pdf("MDSplot.pdf")
plotMDS(y,cex=0.8,col=pal)
dev.off()

# used for Camila's MSD since there were only two groups/one comparison
#pdf("MDSplot.pdf")
#plotMDS(y,cex=0.8,col=c(rep("red",2), rep("blue", 3), rep("red", 2), "blue", "red", "blue", "blue", "red"))
#dev.off()

d<-cpm(y)
write.table(d,"CPM.txt",quote=FALSE,sep="\t",row.names=TRUE)

var_genes <- apply(d, 1, var)
var_50 <- names(sort(var_genes, decreasing=TRUE))[1:50]
cpm_50 <- d[var_50,]

pdf("top50genes.edgeR.heatmap.pdf", height=10)
pheatmap(cpm_50,
         scale="row",
         color=inferno(30),
         fontsize_row = 8,
         main=paste0("Normalized Counts for Top 50 Variable Genes"),
)
dev.off()

designlist <- sapply(grouplist,
                     function(x) {
                       model.matrix(~0+x)
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

ylist <- lapply(designlist,
                function (x) {
                  x <- y
                  })

ylist <- sapply(names(ylist),
                function(x) {
                  ylist[[x]] <- estimateGLMCommonDisp(ylist[[x]], designlist[[x]])
                  ylist[[x]] <- estimateGLMTrendedDisp(ylist[[x]], designlist[[x]])
                  ylist[[x]] <- estimateGLMTagwiseDisp(ylist[[x]], designlist[[x]])
                },
                simplify = FALSE,
                USE.NAMES = TRUE)

qfitlist <- sapply(names(ylist),
                   function(x) {
                     ylist[[x]] <- glmQLFit(ylist[[x]], designlist[[x]])
                   },
                   simplify = FALSE,
                   USE.NAMES = TRUE)

gene.constant <- dim(y)[1]

#### EDGER: DIFFERENTIAL EXPRESSION ####
qlflist <- sapply(names(qfitlist),
                   function (x) {
                     glmQLFTest(contrast = compList[[x]][[2]], glmfit = qfitlist[[x]])},
                  simplify = FALSE,
                  USE.NAMES = TRUE)


top.qlflist <- sapply(names(qlflist),
                      function (x) {
                        topTags(qlflist[[x]], n = gene.constant, p.value = 0.1)},
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#write DEG information for all genes
sapply(names(qlflist),
       function (x) {
         write.table(qlflist[[x]], file = paste0(x, ".all.edgeR.txt"), quote=FALSE, sep = "\t")
       })

#make standard MA plots (MA plots)
sapply(names(qlflist),
       function (x) {
         pdf(file = paste0(x,".edgeR.MAplot.pdf"))
         plotMD(qlflist[[x]], main=x)
         abline(h=c(-1,1), col="blue")
         dev.off()
         },
       simplify = FALSE,
       USE.NAMES = TRUE)

deg.list <- lapply(top.qlflist,
                   as.data.frame)

deg.sig.list <- lapply(deg.list,
                       function(x) x[x$FDR < 0.05,] )

deg.up.list <- lapply(deg.sig.list,
                      function(x) {
                        if (dim(x)[[1]] > 0) {
                          x <- x[x$logFC > 0,]
                          x <- x[order(-x$logFC),]
                        }
                      })

sapply(names(deg.up.list),
       function (x) {
         write.table(deg.up.list[[x]], file = paste0(x, ".deg.up.sig.edgeR.txt"), quote=FALSE, sep='\t')
         })

deg.down.list <- lapply(deg.sig.list,
                      function(x) {
                        if (dim(x)[[1]] > 0) {
                          x <- x[x$logFC < 0,]
                          x <- x[order(-x$logFC),]
                        }
                      })

sapply(names(deg.down.list),
       function (x) {
         write.table(deg.down.list[[x]], file = paste0(x, ".deg.down.sig.edgeR.txt"), quote=FALSE, sep='\t')
       })

deg.all.list <- lapply(qlflist,
                       as.data.frame)
sapply(names(deg.all.list),
       function (x) {
         write.table(deg.list[[x]], file = paste0(x, ".deg.allGenes.edgeR.txt"), quote=FALSE, sep='\t')
       })

#### EXCLUDE MITOCHONDRIAL GENES ####

gcnts$genes <- row.names(gcnts)
gcnts <- gcnts %>%
  filter(!grepl("\\|MT", genes))

gcnts <- subset(gcnts, select=-c(genes))

#### EDGER: MODEL AND MDS PLOT WITHOUT MT GENES ####
y<-DGEList(counts=gcnts)
y<-calcNormFactors(y)

num <- length(gcnts)

pal <- viridis_pal(option = "H", begin=0, end=0.9)(num)

pdf("MDSplot.noMT.pdf")
plotMDS(y,cex=0.8,col=pal)
dev.off()

# used for Camila's MSD since there were only two groups/one comparison
#pdf("MDSplot.noMT.pdf")
#plotMDS(y,cex=0.8,col=c(rep("red",2), rep("blue", 3), rep("red", 2), "blue", "red", "blue", "blue", "red"))
#dev.off()

d<-cpm(y)
write.table(d,"CPM.noMT.txt",quote=FALSE,sep="\t",row.names=TRUE)

var_genes <- apply(d, 1, var)
var_50 <- names(sort(var_genes, decreasing=TRUE))[1:50]
cpm_50 <- d[var_50,]

pdf("top50genes.noMT.edgeR.heatmap.pdf", height=10)
pheatmap(cpm_50,
         scale="row",
         color=inferno(30),
         fontsize_row = 8,
         main=paste0("Normalized Counts for 50 Most Variable Non-Mitochondrial Genes"),
)
dev.off()

designlist <- sapply(grouplist,
                     function(x) {
                       model.matrix(~0+x)
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

ylist <- lapply(designlist,
                function (x) {
                  x <- y
                })

ylist <- sapply(names(ylist),
                function(x) {
                  ylist[[x]] <- estimateGLMCommonDisp(ylist[[x]], designlist[[x]])
                  ylist[[x]] <- estimateGLMTrendedDisp(ylist[[x]], designlist[[x]])
                  ylist[[x]] <- estimateGLMTagwiseDisp(ylist[[x]], designlist[[x]])
                },
                simplify = FALSE,
                USE.NAMES = TRUE)

qfitlist <- sapply(names(ylist),
                   function(x) {
                     ylist[[x]] <- glmQLFit(ylist[[x]], designlist[[x]])
                   },
                   simplify = FALSE,
                   USE.NAMES = TRUE)

gene.constant <- dim(y)[1]

#### EDGER: DIFFERENTIAL EXPRESSION ####
qlflist <- sapply(names(qfitlist),
                  function (x) {
                    glmQLFTest(contrast = compList[[x]][[2]], glmfit = qfitlist[[x]])},
                  simplify = FALSE,
                  USE.NAMES = TRUE)


top.qlflist <- sapply(names(qlflist),
                      function (x) {
                        topTags(qlflist[[x]], n = gene.constant, p.value = 0.1)},
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#write DEG information for all genes
sapply(names(qlflist),
       function (x) {
         write.table(qlflist[[x]], file = paste0(x, ".noMT.all.edgeR.txt"), quote=FALSE, sep = "\t")
       })

#make standard MA plots (MA plots)
sapply(names(qlflist),
       function (x) {
         pdf(file = paste0(x,".noMT.edgeR.MAplot.pdf"))
         plotMD(qlflist[[x]], main=x)
         abline(h=c(-1,1), col="blue")
         dev.off()
       },
       simplify = FALSE,
       USE.NAMES = TRUE)

deg.list <- lapply(top.qlflist,
                   as.data.frame)

deg.sig.list <- lapply(deg.list,
                       function(x) x[x$FDR < 0.05,] )

deg.up.list <- lapply(deg.sig.list,
                      function(x) {
                        if (dim(x)[[1]] > 0) {
                          x <- x[x$logFC > 0,]
                          x <- x[order(-x$logFC),]
                        }
                      })

sapply(names(deg.up.list),
       function (x) {
         write.table(deg.up.list[[x]], file = paste0(x, ".noMT.deg.up.sig.edgeR.txt"), quote=FALSE, sep='\t')
       })

deg.down.list <- lapply(deg.sig.list,
                        function(x) {
                          if (dim(x)[[1]] > 0) {
                            x <- x[x$logFC < 0,]
                            x <- x[order(-x$logFC),]
                          }
                        })

sapply(names(deg.down.list),
       function (x) {
         write.table(deg.down.list[[x]], file = paste0(x, ".noMT.deg.down.sig.edgeR.txt"), quote=FALSE, sep='\t')
       })

deg.all.list <- lapply(qlflist,
                       as.data.frame)
sapply(names(deg.all.list),
       function (x) {
         write.table(deg.list[[x]], file = paste0(x, ".noMT.deg.allGenes.edgeR.txt"), quote=FALSE, sep='\t')
       })
