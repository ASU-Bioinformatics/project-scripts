#### LOAD LIBRARIES ####
pkg1 <- c('NOISeq', 'viridis', "tidyverse", "optparse", "dplyr")
lapply(pkg1,
       function (x) library(x, character.only = TRUE))


#### READ IN ARGUMENTS ####
option_list <- list(
  make_option(c("-d", "--directory"), type="character",
              default="/Volumes/Gencore/analysis_projects/6196658_Sudhindra/comparative_transcriptomics/capensis-only/alignment-short/stringtie_out",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character",
              default="gene.count.matrix.tsv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="comparisons.csv",
              help="comparisons for differential expression")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$directory)

gcnts <-read.delim(opts$genes, header=TRUE, row.names = 1)
gcnts %>%
  select(sort(names(.)))

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = "Name")
comparisons <- comparisons[ order(row.names(comparisons)), ]

colnames(gcnts) <- rownames(comparisons)
gcnts <- head(gcnts, -5)

#define comparisons for NOIseq
compList <- sapply(names(comparisons),
                   function (x) {
                     x <- list(x, c(unique(na.omit(comparisons[,grep(x, colnames(comparisons))]))))
                   },
                   simplify = FALSE,
                   USE.NAMES = TRUE)

genes.filt = lapply(comparisons,
                    function (x) {
                      filtered.data(gcnts, factor=x,
                                    norm = FALSE, method = 1,
                                    cv.cutoff = 100, cpm = 1,
                                    p.adj = "fdr")
                    })

rData.list <- sapply(names(genes.filt),
                     function (x) {
                       NOISeq::readData(data = genes.filt[[x]],
                                        factors = as.data.frame(comparisons[[x]]))
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

rTMM.list <- lapply(rData.list,
                    function (x) {
                      tmm(assayData(x)$exprs, long = 1000, lc = 0)
                    })

rTMM.filt.list <- sapply(names(rTMM.list),
                         function (x) {
                           filtered.data(rTMM.list[[x]], factor=comparisons[[x]],
                                         norm = FALSE, method = 1, cv.cutoff = 95,
                                         cpm = 3, p.adj = "bonferroni")
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

tmmData.list <- sapply(names(rTMM.filt.list),
                       function (x) {
                         NOISeq::readData(data = rTMM.filt.list[[x]],
                                          factors = as.data.frame(comparisons[[x]]))
                       },
                       simplify = FALSE,
                       USE.NAMES = TRUE)

tmmPCA.list <- lapply(tmmData.list,
                      dat,
                      type = "PCA")

sapply(names(tmmPCA.list),
       function (x) {
         pdf(paste0("PCAplot.", x, ".noiseq.pdf"))
         y <- x[1]
         y <- str_split(y, '_')
         y <- as.data.frame(y)
         contrast <- y[[1]][1]
         base <- y[[1]][3]
         xp <- explo.plot(tmmPCA.list[[x]]) +
           scale_color_hue(labels = c(base,contrast))
         print(xp)
         dev.off()
       })

results.noiseq.list <- sapply(names(comparisons),
                              function (x) {
                                if ((table(comparisons[x])[[1]] > 1) && (table(comparisons[x])[[2]] > 1)) {
                                  noiseqbio(tmmData.list[[x]], k = 0.5, norm = "tmm",
                                            factor = "comparisons[[x]]",
                                            conditions = c(-1, 1),
                                            a0per = 0.9, r = 20, lc = 0)
                                } else {
                                  noiseq(tmmData.list[[x]], k=0.5, norm = "tmm", replicates = 'no',
                                         factor = "comparisons[[x]]",
                                         conditions = c(-1, 1))
                                }
                              },
                              simplify = FALSE,
                              USE.NAMES = TRUE)

save(results.noiseq.list, file = "NOIseq_results.rData")

degenes.lists <- sapply(names(results.noiseq.list),
                        function (x) {
                          lists <- list("all" = degenes(results.noiseq.list[[x]], q=0.99, M=NULL),
                                        "up.sig" = degenes(results.noiseq.list[[x]], q=0.99, M="up"),
                                        "down.sig" = degenes(results.noiseq.list[[x]], q=0.99, M="down"))
                          return(lists)
                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE)

sapply(names(degenes.lists),
       function (x) {
         sapply(names(degenes.lists[[x]]),
                function (y) {
                  write.table(degenes.lists[[x]][[y]],
                              file = paste0(x,".deg.",y,".noiseq.txt"),
                              sep = "\t", quote = FALSE)
                })
       })

sapply(names(results.noiseq.list),
       function (x) {
         pdf(file = paste0(x,".noiseq.MDplot.pdf"))
         pm <- DE.plot(results.noiseq.list[[x]], q = 0.99, graphic = "MD", log.scale = TRUE)
         print(pm)
         dev.off()
         pdf(file = paste0(x,".noiseq.explot.pdf"))
         pe <- DE.plot(results.noiseq.list[[x]], q = 0.99, graphic = "expr", log.scale = TRUE)
         print(pe)
         dev.off()
       })
