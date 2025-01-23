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
  make_option(c("-d", "--directory"), type="character", default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_5/short-myxoma-DEGs-2/",
              help="path to count matrix and template information"),
  make_option(c("-g", "--genes"), type="character", default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_5/short-stringtie/gene_count_matrix_myxoma.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character", default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_5/comparisons.csv",
              help="comparisons for differential expression"),
  make_option(c("-l", "--logFCcutoff"), type="integer", default = 1,
              help = "log2 fold change cutoff for stringent tables and plots"),
  make_option(c("-n", "--toolNumber"), type="integer", default = 2,
              help = "minimum number of tools in which a gene is identified as significant DEG"),
  make_option(c("-f", "--functionToGene"), type="character", default = "/Volumes/Gencore/databases/reference_genomes/human/human-myxoma-dual-STAR2.7.10/myxoma.goterms.txt",
              help = "text file with two columns correlating GO terms to genes")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$directory)

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = "Name")

#### READ IN DIFFERENTIAL GENE LISTS ####

deg.up.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".up.0.0.stats.csv"), sep=",", row.names = "gene_id")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".down.0.0.stats.csv"), sep=",", row.names = "gene_id")
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
pkg2 <- c('igraph', 'clusterProfiler', 'enrichplot', 'ggnewscale')
lapply(pkg2,
       function(x) library(x, character.only = TRUE))

library(GO.db)
GO <- as.list(GOTERM)
GO$`GO:0000001`@Term

GO <- as.list(Term(GOTERM))


term.to.genes <- read.delim(opts$functionToGene, header=FALSE, row.names = NULL)

#enricher(
 # gene = deg.down.list$selinexorM159_v_M159$ids,
#  minGSSize = 3,
#  pvalueCutoff = 0.5,
#  pAdjustMethod = "BH",
#  qvalueCutoff = 0.5,
#  TERM2GENE = term.to.genes
#  )

up.cp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           enricher(
                             gene = deg.up.list[[x]]$ids,
                             minGSSize = 3,
                             pvalueCutoff = 0.5,
                             pAdjustMethod = "BH",
                             qvalueCutoff = 1,
                             TERM2GENE = term.to.genes
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
                             qvalueCutoff = 1,
                             TERM2GENE = term.to.genes
                           )
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

up.cp.deg.list <- sapply(names(deg.up.list),
                           function(x) {
                             if (!is.null(up.cp.deg.list[[x]])) {
                               up.cp.deg.list[[x]]@result$Description <- unname(t(as.data.frame(GO[c(up.cp.deg.list[[x]]@result$ID)])))
                               return(up.cp.deg.list[[x]])
                             }
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             if (!is.null(down.cp.deg.list[[x]])) {
                               down.cp.deg.list[[x]]@result$Description <- unname(t(as.data.frame(GO[c(down.cp.deg.list[[x]]@result$ID)])))
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
