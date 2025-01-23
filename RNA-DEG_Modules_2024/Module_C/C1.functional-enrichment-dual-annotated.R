#### READ IN ARGUMENTS ####
library("optparse", "dplyr", "tidyverse")

option_list <- list(
  make_option(c("-d", "--directory"), type="character", default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_5/short-DEGs/deglists",
              help="path to count matrix and template information"),
  make_option(c("-p", "--scripts"), type="character", default="/Volumes/Gencore/shared_scripts/RNAseq",
              help="path to scripts directory containing meta information file"),
  make_option(c("-g", "--genes"), type="character", default="gene_count_matrix_idsonly.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character", default="/Volumes/Gencore/analysis_projects/5663018_Masmudur_RNA/project_5/comparisons.csv",
              help="comparisons for differential expression"),
  make_option(c("-s", "--species"), type="character", default="homo_sapiens",
              help="species of origin"),
  make_option(c("-l", "--logFCcutoff"), type="integer", default = 1,
              help = "log2 fold change cutoff for stringent tables and plots"),
  make_option(c("-n", "--toolNumber"), type="integer", default = 2,
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

deg.up.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".up.0.0.stats.csv"), sep=",")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.data <- sapply(colnames(comparisons),
                      function(x) {
                        x <- read.delim(paste0("merged_deg.", x, ".down.0.0.stats.csv"), sep=",")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#compute the average fold change value across the three tools
#keep only the genes identified by at least the number of specified tools in the options
total.tools <- (opts$toolNumber - ((length(deg.up.data[[1]]) - 3) /2)) * 2

deg.up.list <- sapply(deg.up.data,
                      function(x) {
                        x <- x[rowSums(is.na(x)) <= total.tools, ]
                        x <- x[order(-x$average.padj),]
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.list <- sapply(deg.down.data,
                      function(x) {
                        x <- x[rowSums(is.na(x)) <= total.tools, ]
                        x <- x[order(-x$average.padj),]
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

#### GET GO TERMS FOR DEG GENE LISTS ####
library(biomaRt)
mart <- useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
# human 'universe' or background GO term/gene list
human.results <- getBM(attributes = c("ensembl_gene_id", "go_id", "name_1006", "definition_1006", "namespace_1003"), mart = mart)

library(GO.db, dplyr)
all.go.terms <- dplyr::distinct(toTable(GOTERM)[ , 2:5])
myxoma.terms <- read.delim("myxoma.goterms.txt", sep = '\t', header = FALSE) # this is a reference file I created using myxoma gff file
# myxoma 'universe' or background GO term/gene list
myxoma.results <- merge(x = myxoma.terms, y = all.go.terms, by.x = "V1", by.y = "go_id")

library("clusterProfiler")

enricher(deg.up.list$)

#### FUNCTIONAL ENRICHMENT ANALYSIS: GPROFILER2 ####
library("gprofiler2")
cb.gp.deg.list <- sapply(names(deg.up.list),
                         function(x) {
                           gost(list("increased" = row.names(deg.up.list[[x]]),
                                     "decreased" = row.names(deg.down.list[[x]])),
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

#cb.gp.deg.list = cb.gp.deg.list[-which(sapply(cb.gp.deg.list, is.null))]
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
                           enrichGO(gene = row.names(deg.up.list[[x]]),
                                    OrgDb = species[3],
                                    ont = "ALL",
                                    keyType = 'ENSEMBL',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = 0.05,
                                    qvalueCutoff = 0.05)
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

#save(up.cp.deg.list, file = "clusterProfiler_up_results.rData")
#load(file = "clusterProfiler_up_results.rData")

down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             enrichGO(gene = row.names(deg.down.list[[x]]),
                                      OrgDb = species[3],
                                      ont = "ALL",
                                      keyType = 'ENSEMBL',
                                      pAdjustMethod = "BH",
                                      pvalueCutoff = 0.05,
                                      qvalueCutoff = 0.05)
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

#save(down.cp.deg.list, file = "clusterProfiler_down_results.rData")
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

#### CLUSTERPROFILER KEGG Analysis ####
library(biomaRt)
mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL",
                         dataset = "hsapiens_gene_ensembl",
                         host = "https://www.ensembl.org")

genes.down <- sapply(names(deg.down.list),
                function(x) {
                  getBM(filters = "ensembl_gene_id",
                        attributes = c("ensembl_gene_id","entrezgene_id","uniprot_gn_id"),
                        values = rownames(deg.down.list[[x]]), 
                        mart = mart)
                },
                simplify = FALSE,
                USE.NAMES = TRUE
)

genes.up <- sapply(names(deg.up.list),
                     function(x) {
                       getBM(filters = "ensembl_gene_id",
                             attributes = c("ensembl_gene_id","entrezgene_id","uniprot_gn_id"),
                             values = rownames(deg.up.list[[x]]), 
                             mart = mart)
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE
)

kegg.down.cp.deg.list <- sapply(names(deg.down.list),
                           function(x) {
                             enrichKEGG(gene = genes.down[[x]]$entrezgene_id,
                                        organism = 'hsa',
                                        keyType = 'kegg',
                                        pvalueCutoff = 0.05)
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

kegg.up.cp.deg.list <- sapply(names(deg.up.list),
                                function(x) {
                                  enrichKEGG(gene = genes.up[[x]]$entrezgene_id,
                                             organism = 'hsa',
                                             keyType = 'kegg',
                                             pvalueCutoff = 0.05)
                                },
                                simplify = FALSE,
                                USE.NAMES = TRUE)

lapply(names(kegg.up.cp.deg.list),
       function (y) {
         if (!is.null(kegg.up.cp.deg.list[[y]])) {
           if (dim(kegg.up.cp.deg.list[[y]]@result)[[1]] != 0) {
             plot.bar <- mutate(kegg.up.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore")
             pdf(file = paste0("kegg.cp.up.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(kegg.up.cp.deg.list[[y]], showCategory=20)
             pdf(file = paste0("kegg.cp.up.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)

lapply(names(kegg.down.cp.deg.list),
       function (y) {
         if (!is.null(kegg.down.cp.deg.list[[y]])) {
           if (dim(kegg.down.cp.deg.list[[y]]@result)[[1]] != 0) {
             plot.bar <- mutate(kegg.down.cp.deg.list[[y]], qscore = -log(p.adjust, base=10)) %>% 
               barplot(x="qscore")
             pdf(file = paste0("kegg.cp.down.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(kegg.down.cp.deg.list[[y]], showCategory=20)
             pdf(file = paste0("kegg.cp.down.dotplot.",y,".pdf"))
             print(plot.dot)
             dev.off()
           }
         }
       }
)

kegg.up.cp.deg.list <- sapply(names(kegg.up.cp.deg.list),
                         function (y) {
                           if (!is.null(kegg.up.cp.deg.list[[y]])) {
                             if (dim(kegg.up.cp.deg.list[[y]]@result)[[1]] > 0) {
                               pairwise_termsim(kegg.up.cp.deg.list[[y]])
                             }
                           }
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

kegg.down.cp.deg.list <- sapply(names(kegg.down.cp.deg.list),
                           function (y) {
                             if (!is.null(kegg.down.cp.deg.list[[y]])) {
                               if (dim(kegg.down.cp.deg.list[[y]]@result)[[1]] > 0) {
                                 pairwise_termsim(kegg.down.cp.deg.list[[y]])
                               }
                             }
                           },
                           simplify = FALSE,
                           USE.NAMES = TRUE)

lapply(names(kegg.down.cp.deg.list),
       function (y) {
         if (!is.null(kegg.down.cp.deg.list[[y]])) {
           if (dim(kegg.down.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(kegg.down.cp.deg.list[[y]])
             pdf(file = paste0("kegg.cp.down.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(kegg.up.cp.deg.list),
       function (y) {
         if (!is.null(kegg.up.cp.deg.list[[y]])) {
           if (dim(kegg.up.cp.deg.list[[y]]@result)[[1]] > 1) {
             plot.emap <- emapplot(kegg.up.cp.deg.list[[y]])
             pdf(file = paste0("kegg.cp.up.emapplot.",y,".pdf"))
             print(plot.emap)
             dev.off()
           }
         }
       })

lapply(names(kegg.up.cp.deg.list),
       function (x) {
         if (!is.null(kegg.up.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(kegg.up.cp.deg.list[[x]],
                                      M = deg.up.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(kegg.up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(kegg.up.cp.deg.list[[x]],
                                     M = deg.up.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("kegg.cp.up.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

lapply(names(kegg.down.cp.deg.list),
       function (x) {
         if (!is.null(kegg.down.cp.deg.list[[x]])) {
           plot.cnet.gene <- cnetplot(kegg.down.cp.deg.list[[x]],
                                      M = deg.down.list[[x]],
                                      cex_label_gene = 0.6,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(kegg.down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.6,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(kegg.down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_gene = 0.5,
                                     cex_label_category = 0.7,
                                     node_label = "all")
           pdf(file = paste0("kegg.cp.down.cnetplot.",x,".pdf"))
           print(plot.cnet.gene)
           print(plot.cnet.cat)
           print(plot.cnet.all)
           dev.off()
         }
       })

sapply(names(kegg.up.cp.deg.list),
       function (x) {
         if (!is.null(kegg.up.cp.deg.list[[x]])) {
           write.table(kegg.up.cp.deg.list[[x]]@result, file = paste0("kegg.", x, ".clusterProfiler.up.sig.csv"), quote=FALSE, sep=',')
         }})

sapply(names(kegg.down.cp.deg.list),
       function (x) {
         if (!is.null(kegg.down.cp.deg.list[[x]])) {
           write.table(kegg.down.cp.deg.list[[x]]@result, file = paste0("kegg.", x, ".clusterProfiler.down.sig.csv"), quote=FALSE, sep=',')
         }})


