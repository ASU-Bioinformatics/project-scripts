#### README ####
# This script runs functional enrichment for an organism that exists in the KEGG database

# Inputs include a gene count matrix, differential expression files (all from edgeR, merged significant across all tools),
# pairwise comparison metadata files, and a species key for KEGG organism identification

# The first step reads in all the differential expression files, extracts the IDS from the row names,
# filters by the tool stringency desired (eg, significant according to at least two tools), and sorts by the p-adj
# There will be two lists at this point, for the up-regulated and down-regulated genes respectively

# The second step uses gprofiler2 to create a Manhattan plot of significant GO terms and a corresponding text file

# The third step uses clusterProfiler to run GO functional analysis
# This outputs several types of graphs as well as tabular output

# The fourth step also uses clusterProfiler, but runs KEGG pathway analysis instead
# This outputs gsea plots, pathview figures, and tables for each identified differential KEGG pathway for each pairwise comparison


#### READ IN ARGUMENTS ####
library("optparse", "dplyr", "tidyverse", 'stringr')

option_list <- list(
  make_option(c("-d", "--directory"), type="character", 
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/DEGS-lowqual",
              help="path to merged files and output"),
  make_option(c("-p", "--scripts"), type="character", default="/Volumes/Gencore/shared_scripts/RNAseq",
              help="path to scripts directory containing meta information file"),
  make_option(c("-g", "--genes"), type="character",
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/stringtie-lowqual/gene_count_matrix.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="/Volumes/Gencore/analysis_projects/7232055_Nolz/50M-reads-deadline/comparisons.csv",
              help="comparisons for differential expression"),
  make_option(c("-s", "--species"), type="character", default="homo_sapiens",
              help="species of origin"),
  make_option(c("-l", "--logFCcutoff"), type="integer", default = 0,
              help = "log2 fold change cutoff for stringent tables and plots"),
  make_option(c("-n", "--toolNumber"), type="integer", default = 1,
              help = "minimum number of tools in which a gene is identified as significant DEG")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$scripts)

spc <- opts$species
species.meta <- read.csv("species-meta.csv", header=TRUE, row.names = "Data")
species <- species.meta[[spc]]

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

#### FUNCTIONAL ENRICHMENT ANALYSIS: GPROFILER2 ####
library("gprofiler2")
cb.gp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           gost(list("increased" = deg.up.list[[x]]$ids,
                                     "decreased" = deg.down.list[[x]]$ids),
                                organism = species[1],
                                significant = TRUE,
                                ordered_query = TRUE,
                                sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG"))
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

save(cb.gp.deg.list, file = "gprofiler_results.rData")


cb.gp.deg.list <- lapply(cb.gp.deg.list,
                         function(x) {
                           sorted <- x$result[order(x$result$p_value),]
                           return(list(result = sorted,
                                       meta = x$meta))})

cb.gp.deg.list = cb.gp.deg.list[-which(sapply(cb.gp.deg.list, is.null))]
lapply(names(cb.gp.deg.list),                  
       function(x) {
         publish_gostplot(gostplot(cb.gp.deg.list[[x]], interactive = FALSE),
                          highlight_terms = cb.gp.deg.list[[x]]$result$term_id[1:5],
                          width = 12,
                          filename = paste0("gostplot-", x, ".png"))
       })

cb.gp.results <- sapply(names(cb.gp.deg.list),
                        function (x) as.data.frame(cb.gp.deg.list[[x]]$result[1:13]),
                        simplify = FALSE,
                        USE.NAMES = TRUE)

sapply(names(cb.gp.results),
       function (x) write.table(cb.gp.results[[x]], file = paste0(x, ".gprofiler2.sig.txt"), quote=FALSE))

#### FUNCTIONAL ENRICHMENT ANALYSIS: CLUSTERPROFILER ####
detach("package:edgeR", unload=TRUE)
pkg2 <- c('igraph', 'clusterProfiler', 'enrichplot', 'ggnewscale', species[[3]])
lapply(pkg2,
       function(x) library(x, character.only = TRUE))

up.cp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           enrichGO(gene = deg.up.list[[x]]$ids,
                                    OrgDb = species[3],
                                    ont = "ALL",
                                    keyType = 'ENSEMBL',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

save(up.cp.deg.list, file = "clusterProfiler_up_results.rData")
#load(file = "clusterProfiler_up_results.rData")

down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             enrichGO(gene = deg.down.list[[x]]$ids,
                                      OrgDb = species[3],
                                      ont = "ALL",
                                      keyType = 'ENSEMBL',
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05)
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

save(down.cp.deg.list, file = "clusterProfiler_down_results.rData")
#load(file = "clusterProfiler_down_results.rData")

lapply(names(up.cp.deg.list),
       function (y) {
         if (!is.null(up.cp.deg.list[[y]])) {
          if (dim(up.cp.deg.list[[y]]@result)[[1]] != 0) {
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
           if (dim(down.cp.deg.list[[y]]@result)[[1]] != 0) {
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

rd.up.cp.deg.list <- sapply(names(up.cp.deg.list),
                            function (x) {
                              if (!is.null(up.cp.deg.list[[x]])) {
                                setReadable(up.cp.deg.list[[x]],
                                            OrgDb = species[3],
                                            keyType = "ENSEMBL")
                              }},
                            simplify = FALSE,
                            USE.NAMES = TRUE)

rd.down.cp.deg.list <- sapply(names(down.cp.deg.list),
                              function (x) {
                                if (!is.null(down.cp.deg.list[[x]])) {
                                  setReadable(down.cp.deg.list[[x]],
                                              OrgDb = species[3],
                                              keyType = "ENSEMBL")
                                }},
                              simplify = FALSE,
                              USE.NAMES = TRUE)

lapply(names(rd.up.cp.deg.list),
       function (x) {
         if (!is.null(rd.up.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(rd.up.cp.deg.list[[x]],
                                      M = deg.up.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(rd.up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(rd.up.cp.deg.list[[x]],
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

lapply(names(rd.down.cp.deg.list),
       function (x) {
         if (!is.null(rd.down.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(rd.down.cp.deg.list[[x]],
                                      M = deg.down.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(rd.down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(rd.down.cp.deg.list[[x]],
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
         if (!is.null(rd.up.cp.deg.list[[x]])) {
           write.table(up.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.up.sig.csv"), quote=FALSE, sep = ',')
         }})

sapply(names(down.cp.deg.list),
       function (x) {
         if (!is.null(rd.down.cp.deg.list[[x]])) {
           write.table(down.cp.deg.list[[x]]@result, file = paste0("go.", x, ".clusterProfiler.down.sig.csv"), quote=FALSE, sep = ',')
         }})

#### CLUSTERPROFILER KEGG ANALYSIS ####
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org")

all.genes <- sapply(names(comparisons),
                    function(x) {
                      x <- read.delim(paste0(x, '.all.edgeR.txt'))
                      x$ids <- as.data.frame(stringr::str_split_fixed(rownames(x), stringr::fixed("|"), 2))$V1
                      return(x)
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE)

all.genes.entrez <- sapply(names(all.genes),
                           function(x) {
                             getBM(filters = "ensembl_gene_id",
                                   attributes = c("ensembl_gene_id","entrezgene_id"),
                                   values = all.genes[[x]]$ids, 
                                   mart = mart)
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

all.genes.entrez <- sapply(names(all.genes),
                           function(x) {
                             merge(all.genes[[x]], all.genes.entrez[[x]],
                                   by.x = 'ids', by.y = 'ensembl_gene_id')
                           }, 
                           simplify = FALSE,
                           USE.NAMES = TRUE)

all.genes.entrez <- lapply(all.genes.entrez,
                           function(x) {
                             x <- x[which(duplicated(x$entrezgene_id) == F), ]
                           })

foldchanges <- sapply(names(all.genes.entrez),
                      function(x) {
                        print(x)
                        y <- all.genes.entrez[[x]]$logFC
                        names(y) <- all.genes.entrez[[x]]$entrezgene_id
                        y <- sort(y, decreasing = TRUE)
                        return(y)
                      }, 
                      simplify = FALSE,
                      USE.NAMES = TRUE)

gseaKEGGs <- sapply(names(foldchanges),
                    function(x) {
                      gseaKEGG <- gseKEGG(geneList = foldchanges[[x]], # ordered named vector of fold changes (Entrez IDs are the associated names)
                                          organism = "hsa", # supported organisms listed below
                                          nPermSimple = 10000, # default number permutations
                                          minGSSize = 20, # minimum gene set size (# genes in set) - change to test more sets or recover sets with fewer # genes
                                          pvalueCutoff = 0.05, # padj cutoff value
                                          verbose = FALSE)
                    },
                    simplify=FALSE,
                    USE.NAMES=TRUE)

gseaKEGGresults <- sapply(names(gseaKEGGs),
                          function(x) {
                            x <- gseaKEGGs[[x]]@result
                          },
                          simplify = FALSE,
                          USE.NAMES = TRUE)

sapply(names(gseaKEGGs),
       function (x) {
         if (!is.null(gseaKEGGresults[[x]])) {
           write.table(gseaKEGGresults[[x]], file = paste0("kegg.", x, ".clusterProfiler.sig.csv"), quote=FALSE, sep = ',')
         }})

sapply(names(gseaKEGGs),
       function(x) {
         pdf(file = paste0("gsea-plot.", x, ".pdf"))
          for (id in gseaKEGGs[[x]]@result$ID) {
            gs <- gseaplot(gseaKEGGs[[x]],
                          geneSetID = id,
                          title = paste0(id, ": ", gseaKEGGs[[x]]@result[id,2]))
           print(gs)
          }
         dev.off()
       })

detach("package:dplyr", unload=TRUE) # first unload dplyr to avoid conflicts
library('pathview')

## Output images for a single significant KEGG pathway
sapply(names(gseaKEGGs),
       function(x) {
         for (id in gseaKEGGs[[x]]@result$ID) {
           skip_to_next <- FALSE
           tryCatch(
             pathview(gene.data = foldchanges[[x]],
                      pathway.id = id,
                      out.suffix = x,
                      species = "hsa",
                      limit = list(gene = 2, # value gives the max/min limit for foldchanges
                                  cpd = 1)),
             error = function(e) { skip_to_next <<- TRUE})
           if (skip_to_next) { next }
             
         }
       }
)


