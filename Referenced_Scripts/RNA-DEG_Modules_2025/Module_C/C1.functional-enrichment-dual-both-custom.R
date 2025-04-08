#### README ####
# This script runs functional enrichment for an organism that isn't in biomart or kegg (custom functional information)

# Inputs include a gene count matrix, differential expression files (all from edgeR, merged significant across all tools),
# pairwise comparison metadata files, and a custom GO term to gene text file

# The first step reads in all the differential expression files, extracts the IDS from the row names,
# filters by the tool stringency desired (eg, significant according to at least two tools), and sorts by the p-adj
# There will be two lists at this point, for the up-regulated and down-regulated genes respectively

# The third step uses clusterProfiler to run GO functional analysis, using provided custom GO term to gene information
# This outputs several types of graphs as well as tabular output

#### READ IN ARGUMENTS ####
library("optparse", "dplyr", "tidyverse", 'stringr')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default="/Volumes/Gencore/analysis_projects/7763550_Junior_RNA/degs-hare-short/degs_withMT",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character", default="/Volumes/Gencore/analysis_projects/7763550_Junior_RNA/stringtie-hare-short/gene_count_matrix_idsonly.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character", default="/Volumes/Gencore/analysis_projects/7763550_Junior_RNA/comparisons.csv",
              help="comparisons for differential expression"),
  make_option(c("-l", "--logFCcutoff"), type="integer", default = 1,
              help = "log2 fold change cutoff for stringent tables and plots"),
  make_option(c("-n", "--toolNumber"), type="integer", default = 2,
              help = "minimum number of tools in which a gene is identified as significant DEG"),
  make_option(c("-1", "--functionToGene1"), type="character", default = "/Volumes/Gencore/databases/reference_genomes/hare/hare-myxoma-dual-STAR2.7.10a/hare-goterms.txt",
              help = "text file with two columns correlating GO terms to genes"),
  make_option(c("-2", "--functionToGene2"), type="character", default = "/Volumes/Gencore/databases/reference_genomes/hare/hare-myxoma-dual-STAR2.7.10a/myxoma.goterms.txt",
              help = "text file with two columns correlating GO terms to genes")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$directory)

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = "Name")

#### READ IN DIFFERENTIAL GENE LISTS ####

deg.up.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".up.1.0.stats.csv"), sep=",", row.names = "gene_id")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".down.1.0.stats.csv"), sep=",", row.names = "gene_id")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#compute the average fold change value across the three tools
#keep only the genes identified by at least the number of specified tools in the options
total.tools <- ((((length(deg.up.data[[1]]) - 2) /2) - opts$toolNumber) * 2)

deg.up.list <- sapply(deg.up.data,
                      function(x) {
                        x <- x[rowSums(is.na(x)) <= total.tools, ]
                        x <- x[order(-x$average.padj),]
                        x$ids <- as.data.frame(stringr::str_split_fixed(rownames(x), stringr::fixed("|"), 2))$V1
                        return(x)
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.list <- sapply(deg.down.data,
                      function(x) {
                        x <- x[rowSums(is.na(x)) <= total.tools, ]
                        x <- x[order(-x$average.padj),]
                        x$ids <- as.data.frame(stringr::str_split_fixed(rownames(x), stringr::fixed("|"), 2))$V1
                        return(x)
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#### FUNCTIONAL ENRICHMENT ANALYSIS: CLUSTERPROFILER ####
detach("package:edgeR", unload=TRUE)
pkg2 <- c('igraph', 'clusterProfiler', 'enrichplot', 'ggnewscale', 'GO.db')
lapply(pkg2,
       function(x) library(x, character.only = TRUE))

all.go.terms <- dplyr::distinct(toTable(GOTERM)[ , 2:5])
myxoma.terms <- read.delim(opts$functionToGene2, sep = '\t', header = FALSE) # this is a reference file I created using myxoma gff file
# myxoma 'universe' or background GO term/gene list
myxoma.results <- merge(x = myxoma.terms, y = all.go.terms, by.x = "V1", by.y = "go_id")

hare.terms <- read.delim(opts$functionToGene1, sep = '\t', header = FALSE) # this is a reference file I created using hare gaf file
# hare 'universe' or background GO term/gene list
hare.results <- merge(x = hare.terms, y = all.go.terms, by.x = "V1", by.y = "go_id")
hare.go <- hare.results[ , c(1,3:5)]
hare.go <- unique(hare.go)

GO <- as.list(GOTERM)
GO$`GO:0000001`@Term

GO <- as.list(Term(GOTERM))


#### SPECIES 1 ####
up.cp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           enricher(
                             gene = deg.up.list[[x]]$ids,
                             minGSSize = 3,
                             pvalueCutoff = 0.01,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             TERM2GENE = hare.terms
                           )
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(deg.down.list),
                         function(x) {
                           enricher(
                             gene = deg.down.list[[x]]$ids,
                             minGSSize = 3,
                             pvalueCutoff = 0.01,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             TERM2GENE = hare.terms
                           )
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

up.cp.deg.list <- sapply(names(deg.up.list),
                           function(x) {
                            if (!is.null(up.cp.deg.list[[x]])) {
                              up.cp.deg.list[[x]]@result <- merge(up.cp.deg.list[[x]]@result, 
                                                                  hare.go[hare.go$V1 %in% up.cp.deg.list[[x]]@result$ID, 1:2], 
                                                                  by.x = "ID", by.y = "V1", sort=FALSE)
                              up.cp.deg.list[[x]]@result$Description <- up.cp.deg.list[[x]]@result$Term
                               return(up.cp.deg.list[[x]])
                             }
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(deg.down.list),
                         function(x) {
                           if (!is.null(down.cp.deg.list[[x]])) {
                             down.cp.deg.list[[x]]@result <- merge(down.cp.deg.list[[x]]@result, 
                                                                 hare.go[hare.go$V1 %in% down.cp.deg.list[[x]]@result$ID, 1:2], 
                                                                 by.x = "ID", by.y = "V1", sort=FALSE)
                             down.cp.deg.list[[x]]@result$Description <- down.cp.deg.list[[x]]@result$Term
                             return(down.cp.deg.list[[x]])
                           }
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

lapply(names(up.cp.deg.list),
       function (y) {
         if (!is.null(up.cp.deg.list[[y]])) {
          if (dim(up.cp.deg.list[[y]])[[1]] != 0) {
            plot.bar <- mutate(up.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
              barplot(x="qscore", font.size = 8)
            pdf(file = paste0("go.cp.up.barplot.",y,".pdf"))
            print(plot.bar)
            dev.off()
            plot.dot <- dotplot(up.cp.deg.list[[y]], showCategory=15, font.size = 8)
            pdf(file = paste0("go.cp.up.dotplot.",y,".pdf"))
            print(plot.dot)
            dev.off()
            }
         }
         }
       )

lapply(names(down.cp.deg.list),
       function (y) {
         if (!is.null(down.cp.deg.list[[y]])) {
           if (dim(down.cp.deg.list[[y]])[[1]] != 0) {
             plot.bar <- mutate(down.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore", font.size = 8)
             pdf(file = paste0("go.cp.down.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(down.cp.deg.list[[y]], showCategory=15, font.size = 8)
             pdf(file = paste0("go.cp.down.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)





up.cp.deg.list <- sapply(names(up.cp.deg.list),
                         function (y) {
                           if (!is.null(up.cp.deg.list[[y]])) {
                            if (dim(up.cp.deg.list[[y]]@result)[[1]] > 0) {
                              pairwise_termsim(up.cp.deg.list[[y]])
                            }
                          }
                          },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(down.cp.deg.list),
                           function (y) {
                             if (!is.null(down.cp.deg.list[[y]])) {
                              if (dim(down.cp.deg.list[[y]]@result)[[1]] > 0) {
                                pairwise_termsim(down.cp.deg.list[[y]])
                              }
                            }
                            },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

lapply(names(down.cp.deg.list),
       function (y) {
         if (!is.null(down.cp.deg.list[[y]])) {
           if (dim(down.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(down.cp.deg.list[[y]])
             pdf(file = paste0("go.cp.down.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(up.cp.deg.list),
       function (y) {
         if (!is.null(up.cp.deg.list[[y]])) {
           if (dim(up.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(up.cp.deg.list[[y]])
             pdf(file = paste0("go.cp.up.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(up.cp.deg.list),
       function (x) {
         if (!is.null(up.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(up.cp.deg.list[[x]],
                                      M = deg.up.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("go.cp.up.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

lapply(names(down.cp.deg.list),
       function (x) {
         if (!is.null(down.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(down.cp.deg.list[[x]],
                                      M = deg.down.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("go.cp.down.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

sapply(names(up.cp.deg.list),
       function (x) {
         if (!is.null(up.cp.deg.list[[x]])) {
           write.table(up.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.up.sig.tsv"), quote=FALSE, sep = '\t', row.names = FALSE)
         }})

sapply(names(down.cp.deg.list),
       function (x) {
         if (!is.null(down.cp.deg.list[[x]])) {
           write.table(down.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.down.sig.tsv"), quote=FALSE, sep = '\t', row.names = FALSE)
         }})

#### SPECIES 2 ####
up.cp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           enricher(
                             gene = deg.up.list[[x]]$ids,
                             minGSSize = 3,
                             pvalueCutoff = 0.5,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.5,
                             TERM2GENE =  myxoma.terms
                           )
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             enricher(
                               gene = deg.down.list[[x]]$ids,
                               minGSSize = 3,
                               pvalueCutoff = 0.5,
                               pAdjustMethod = "BH",
                               qvalueCutoff = 0.5,
                               TERM2GENE =  myxoma.terms
                             )
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

is.NullOb <- function(x) is.null(x) | all(sapply(x, is.null))
rmNullObs <- function(x) {
  x <- Filter(Negate(is.NullOb), x)
  lapply(x, function(x) if (is.list(x)) rmNullObs(x) else x)
}

test <- GO[c(up.cp.deg.list$Laus_v_Mock.noMT@result$ID)]
test1 <- rmNullObs(test)
test2 <- t(as.data.frame(test1))
test2$GOs <- rownames(test2)

up.cp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           if (!is.null(up.cp.deg.list[[x]])) {
                             up.cp.deg.list[[x]]@result <- merge(up.cp.deg.list[[x]]@result, 
                                                                  myxoma.go[ myxoma.go$V1 %in% up.cp.deg.list[[x]]@result$ID, 1:2], 
                                                                 by.x = "ID", by.y = "V1", sort=FALSE)
                             up.cp.deg.list[[x]]@result$Description <- up.cp.deg.list[[x]]@result$Term
                             return(up.cp.deg.list[[x]])
                           }
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             if (!is.null(down.cp.deg.list[[x]])) {
                               down.cp.deg.list[[x]]@result <- merge(down.cp.deg.list[[x]]@result, 
                                                                      myxoma.go[ myxoma.go$V1 %in% down.cp.deg.list[[x]]@result$ID, 1:2], 
                                                                     by.x = "ID", by.y = "V1", sort=FALSE)
                               down.cp.deg.list[[x]]@result$Description <- down.cp.deg.list[[x]]@result$Term
                               return(down.cp.deg.list[[x]])
                             }
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

lapply(names(up.cp.deg.list),
       function (y) {
         if (!is.null(up.cp.deg.list[[y]])) {
           if (dim(up.cp.deg.list[[y]])[[1]] != 0) {
             plot.bar <- mutate(up.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore")
             pdf(file = paste0("go.cp.up.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(up.cp.deg.list[[y]], showCategory=20)
             pdf(file = paste0("go.cp.up.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)

lapply(names(down.cp.deg.list),
       function (y) {
         if (!is.null(down.cp.deg.list[[y]])) {
           if (dim(down.cp.deg.list[[y]])[[1]] != 0) {
             plot.bar <- mutate(down.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore")
             pdf(file = paste0("go.cp.down.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(down.cp.deg.list[[y]], showCategory=20)
             pdf(file = paste0("go.cp.down.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)





up.cp.deg.list <- sapply(names(up.cp.deg.list),
                         function (y) {
                           if (!is.null(up.cp.deg.list[[y]])) {
                             if (dim(up.cp.deg.list[[y]]@result)[[1]] > 0) {
                               pairwise_termsim(up.cp.deg.list[[y]])
                             }
                           }
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(down.cp.deg.list),
                           function (y) {
                             if (!is.null(down.cp.deg.list[[y]])) {
                               if (dim(down.cp.deg.list[[y]]@result)[[1]] > 0) {
                                 pairwise_termsim(down.cp.deg.list[[y]])
                               }
                             }
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

lapply(names(down.cp.deg.list),
       function (y) {
         if (!is.null(down.cp.deg.list[[y]])) {
           if (dim(down.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(down.cp.deg.list[[y]])
             pdf(file = paste0("go.cp.down.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(up.cp.deg.list),
       function (y) {
         if (!is.null(up.cp.deg.list[[y]])) {
           if (dim(up.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(up.cp.deg.list[[y]])
             pdf(file = paste0("go.cp.up.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(up.cp.deg.list),
       function (x) {
         if (!is.null(up.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(up.cp.deg.list[[x]],
                                      M = deg.up.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("go.cp.up.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

lapply(names(down.cp.deg.list),
       function (x) {
         if (!is.null(down.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(down.cp.deg.list[[x]],
                                      M = deg.down.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("go.cp.down.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

sapply(names(up.cp.deg.list),
       function (x) {
         if (!is.null(up.cp.deg.list[[x]])) {
           write.table(up.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.up.sig.csv"), quote=FALSE, sep = ',', row.names = FALSE)
         }})

sapply(names(down.cp.deg.list),
       function (x) {
         if (!is.null(down.cp.deg.list[[x]])) {
           write.table(down.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.down.sig.csv"), quote=FALSE, sep = ',', row.names = FALSE)
         }})
