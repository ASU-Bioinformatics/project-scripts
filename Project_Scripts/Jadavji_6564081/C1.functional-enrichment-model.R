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
library("dplyr", "tidyverse", 'stringr', "biomaRt")

species.meta <- read.csv("/Volumes/Gencore/shared_scripts/RNAseq/species-meta.csv", header=TRUE, row.names = "Data")
species <- species.meta[["homo_sapiens"]]

setwd("/Volumes/Gencore/analysis_projects/completed_projects/6564081_Jadavji_Spatial/differential-expression")

#### LOAD BIOMART ####

cluster.expression <- read.delim("cortex.integrated.feature.expression.by.cluster.csv", sep=",", row.names = "id")

ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
cluster.ensembls <- getBM(attributes=c('ensembl_gene_id', 'hgnc_symbol', 'entrezgene_id'),
      filters = 'hgnc_symbol',
      values = cluster.expression$gene,
      mart = ensembl)

cluster.found <- merge(cluster.expression, cluster.ensembls, by.x = 'gene', by.y = 'hgnc_symbol')

cluster.found <- unique(cluster.found)

#### READ IN DIFFERENTIAL GENE LISTS ####

names <- sapply(c(0:13), function(x) {paste0('cluster.', x)})
deg.up.data <- sapply(setNames(names, names),
                      function(x) {
                        index <- which(names == x) - 1
                        x <- cluster.found %>%
                         filter( cluster == index, avg_log2FC > 0)
                        
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

deg.down.data <- sapply(setNames(names, names),
                        function(x) {
                          index <- which(names == x) - 1
                          x <- cluster.found %>%
                            filter( cluster == index, avg_log2FC < 0)
                          
                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE)

#### FUNCTIONAL ENRICHMENT ANALYSIS: GPROFILER2 ####
library("gprofiler2")
cb.gp.deg.list <- sapply(names(deg.up.data),
                         function(x) {
                           gost(list("increased" = deg.up.data[[x]]$ensembl_gene_id,
                                     "decreased" = deg.down.data[[x]]$ensembl_gene_id),
                                organism = species[1],
                                significant = TRUE,
                                ordered_query = TRUE,
                                sources = c("GO:MF", "GO:CC", "GO:BP", "KEGG"))
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

cb.gp.deg.list <- lapply(cb.gp.deg.list,
                         function(x) {
                           sorted <- x$result[order(x$result$p_value),]
                           return(list(result = sorted,
                                       meta = x$meta))})

#cb.gp.deg.list = cb.gp.deg.list[-which(sapply(cb.gp.deg.list, is.null))]
lapply(names(cb.gp.deg.list),                  
       function(x) {
         print(x)
         name <- paste0("gostplot-", x, ".png")
         publish_gostplot(gostplot(cb.gp.deg.list[[x]], interactive = FALSE),
                          highlight_terms = cb.gp.deg.list[[x]]$result$term_id[1:5],
                          width = 12, height = 8,
                          filename = name)
       })

cb.gp.results <- sapply(names(cb.gp.deg.list),
                        function (x) as.data.frame(cb.gp.deg.list[[x]]$result[1:13]),
                        simplify = FALSE,
                        USE.NAMES = TRUE)

sapply(names(cb.gp.results),
       function (x) write.table(cb.gp.results[[x]], file = paste0(x, ".gprofiler2.sig.txt"), quote=FALSE, sep="\t"))

#### FUNCTIONAL ENRICHMENT ANALYSIS: CLUSTERPROFILER ####
pkg2 <- c('igraph', 'clusterProfiler', 'enrichplot', 'ggnewscale', species[[3]])
lapply(pkg2,
       function(x) library(x, character.only = TRUE))

up.cp.deg.list <- sapply(names(deg.up.data),
                         function(x) {
                           enrichGO(gene = deg.up.data[[x]]$ensembl_gene_id,
                                    OrgDb = species[3],
                                    ont = "ALL",
                                    keyType = 'ENSEMBL',
                                    pAdjustMethod = "BH",
                                    pvalueCutoff = .05,
                                    qvalueCutoff = .05)
                         },
                         simplify = FALSE,
                         USE.NAMES = TRUE)

save(up.cp.deg.list, file = "clusterProfiler_up_results.rData")
#load(file = "clusterProfiler_up_results.rData")

down.cp.deg.list <- sapply(names(deg.down.data),
                           function(x) {
                             enrichGO(gene = deg.down.data[[x]]$ensembl_gene_id,
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
              barplot(x="qscore", title = paste0("Enriched Functions in ", y, " vs All Other Clusters"))
            pdf(file = paste0("go.cp.up.barplot.",y,".pdf"))
            print(plot.bar)
            dev.off()
            plot.dot <- dotplot(up.cp.deg.list[[y]], showCategory=20,
                                font.size = 6, title = paste0("Enriched Functions in ", y, " vs All Other Clusters"))
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
               barplot(x="qscore", title = paste0("Depleted Functions in ", y, " vs All Other Clusters"))
             pdf(file = paste0("go.cp.down.barplot.",y,".pdf"))
             print(plot.bar)
             dev.off()
             plot.dot <- dotplot(down.cp.deg.list[[y]], showCategory=20,
                                 font.size = 6, title = paste0("Depleted Functions in ", y, " vs All Other Clusters"))
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
             plot.emap <- emapplot(down.cp.deg.list[[y]], cex_line = 0.3, cex_label_category = 0.6,
                                   title = paste0("Depleted Functions in ", y, " vs All Other Clusters"))
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
             plot.emap <- emapplot(up.cp.deg.list[[y]], cex_line = 0.3, cex_label_category = 0.6, repel = TRUE,
                                   title = paste0("Enriched Functions in ", y, " vs All Other Clusters"))
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
                                      M = deg.up.data[[x]],
                                      cex_label_gene = 0.5,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(rd.up.cp.deg.list[[x]],
                                     M = deg.up.data[[x]],
                                     cex_label_category = 0.7,
                                     node_label = "category")
           plot.cnet.all <- cnetplot(rd.up.cp.deg.list[[x]],
                                     M = deg.up.data[[x]],
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
                                      cex_label_gene = 0.5,
                                      node_label = "gene")
           plot.cnet.cat <- cnetplot(rd.down.cp.deg.list[[x]],
                                     M = deg.down.list[[x]],
                                     cex_label_category = 0.7,
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
deg.all.data <- sapply(setNames(names, names),
                      function(x) {
                        index <- which(names == x) - 1
                        x <- cluster.found %>%
                          filter( cluster == index )
                        
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

all.genes.entrez <- lapply(deg.all.data,
                           function(x) {
                             x <- x[which(duplicated(x$entrezgene_id) == F), ]
                           })

foldchanges <- sapply(names(all.genes.entrez),
                      function(x) {
                        print(x)
                        y <- all.genes.entrez[[x]]$avg_log2FC
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


