library("gplots")
library(tidyverse)
library(vegan)
library(ggstatsplot)

setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations")

keggs <- read.delim("ko_map_module_completeness.tab", sep="\t")

kegg.matrix <- keggs[,c(2:33)]
rownames(kegg.matrix) <- kegg.matrix[,1]
kegg.matrix <- kegg.matrix[,2:32]
kegg.matrix$variance <- apply(kegg.matrix[,2:31],1,var)
kegg.matrix <- kegg.matrix[order(kegg.matrix$variance, decreasing = TRUE),]
kegg.modules <- kegg.matrix[1:100,]$pathway.group
kegg.matrix <- as.matrix(kegg.matrix[1:100,2:31])


clustRowBar <- viridis::turbo(length(unique(kegg.modules)), begin=0, end=1)

kegg.modules <- as.data.frame(kegg.modules)
kegg.modules$kegg.modules <- as.factor(kegg.modules$kegg.modules)
kegg.modules$kegg.moduleNames <- as.factor(kegg.modules$kegg.modules)
levels(kegg.modules$kegg.modules) <- as.factor(clustRowBar)
clustRowBar <- as.vector(kegg.modules$kegg.modules)
labels <- c("55", "5A", "4", "4B", "P5", "3C", "P7", "1E", "20",
            "P1", "2D", "P2", "16", "57")

labels <- c("MB-001", "MB-003", "MB-004", "MB-005", "MB-006", "MB-011",
            "MB-012", "MB-014", "MB-016", "MB-017", "MB-018", "MB-020",
            "MB-021", "MB-023", "MB-024", "MB-025", "MB-028", "MB-033",
            "MB-035", "MB-037", "MB-038", "MB-039", "MB-040", "MB-044",
            "MB-045", "MB-047", "MB-049", "MB-050", "MB-053", "MB-055")

pdf("kegg.heatmap.pdf", height=24, width=20)
par(mar = c(2, 2, 16, 2),                                  # Specify par parameters
    xpd = TRUE)
heatmap.2(kegg.matrix[1:100,1:30],
          Rowv=TRUE, 
          Colv=TRUE,
          col=redgreen(100),
          scale="col",
          dendrogram = "col",
          #margins = c(1, 1),
          keysize=0.6,
          cexCol = 1,
          cexRow = 0.9,
          labRow = rownames(kegg.matrix[1:100,1:14]),
          labCol = labels,
          main = "100 Most Variable Kegg Modules by Completeness",
          RowSideColors = clustRowBar,
          trace = "none",
          margins=c(6,44))
#legend(x="top", inset = c(-0.06, -0.12), cex=0.9,
 #      legend = levels(kegg.modules$kegg.moduleNames),
  #     fill = clustRowBar, ncol = 3)
dev.off()

#### PCA PLOT BY GENE COUNT MATRIX ####
metadata <- read.delim("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/metadata.txt", sep="\t")
metadata <- metadata[order(metadata$SampleID),]

setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations")
genecount_df <- read.delim("gene.count.matrix.tsv", sep="\t")
genecount_df <- genecount_df[,order(colnames(genecount_df))]

#colnames(genecount_df) <- c(metadata$Sample.ID)

rownames(metadata) <- metadata$Sample.ID
metadata <- metadata[,-1]

gene_mat = genecount_df |> 
  column_to_rownames("transcript_id") |> 
  as.matrix() |>
  t()

dist_mat = vegdist(gene_mat)

cmd_res = cmdscale(dist_mat, 
                   k = (nrow(gene_mat) - 1),
                   eig = TRUE)

pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

pcoa_meta = bind_cols(pcoa_df, metadata)

#### PLOT PRINCIPAL COMPONENTS COLORED BY METADATA ####
sc <- viridis::scale_color_viridis("viridis", discrete = FALSE)
plots <- lapply(colnames(metadata),
                function (n) {
                  if (is.character(pcoa_meta[[n]])) {
                    m <- pcoa_meta[[n]]
                    p <- ggplot(pcoa_meta,
                                aes(x=PC1, y=PC2, color=m)) +
                      geom_point(size=3) +
                      labs(color = as.character(n), 
                           title = paste0("Principal Coordinates by ", n)) +
                      theme_light() +
                      scale_color_brewer(palette = "Set1")
                    ggsave(filename = paste0("pcoa_plot_", n, ".png"), 
                           plot = p,
                           width = 6,
                           height = 4)}
                  else if (is.integer(pcoa_meta[[n]])) {
                    m <- pcoa_meta[[n]]
                    p <- ggplot(pcoa_meta,
                                aes(x=PC1, y=PC2, color=m)) +
                      geom_point(size=3) +
                      sc +
                      labs(color = as.character(n), 
                           title = paste0("Principal Coordinates by ", n)) +
                      theme_light()
                    ggsave(filename = paste0("pcoa_plot_", n, ".png"), 
                           plot = p,
                           width = 6,
                           height = 4)}
                })


wa_data = wascores(cmd_res$points[,1:2], gene_mat) |>
  as_tibble(rownames = 'geneID')

wa_data

