#### READ IN ARGUMENTS ####
library("optparse")

option_list <- list(
  make_option(c("-d", "--directory"), type="character",
              default="/Volumes/Gencore/analysis_projects/6368449_Hipschman",
              help="path to count matrix and template information"),
  make_option(c("-p", "--scripts"), type="character",
              default="/Volumes/Gencore/shared_scripts/RNAseq",
              help="path to scripts directory containing meta information file"),
  make_option(c("-g", "--genes"), type="character",
              default="gene_count.csv",
              help="count matrix file name"),
  make_option(c("-c", "--comparisons"), type="character",
              default="comparisons.csv",
              help="comparisons for differential expression"),
  make_option(c("-s", "--species"), type="character", default="homo_sapiens",
              help="species of origin"),
  make_option(c("-l", "--logFCcutoff"), type="integer", default = 1,
              help = "log2 fold change cutoff for stringent tables and plots")
)

opt_parser <- OptionParser(option_list = option_list)
opts <- parse_args(opt_parser)

setwd(opts$scripts)

spc <- opts$species
species.meta <- read.csv("species-meta.csv", header=TRUE, row.names = "Data")
species <- species.meta[[spc]]

setwd(opts$directory)

gcnts <-read.csv(opts$genes, header=TRUE, row.names = "gene_id")

comparisons <- read.csv(opts$comparisons, header=TRUE, row.names = "Name")

comparisons[comparisons == -1] <- 2

#ensure that data in the gene count matrix is ordered in the same way as the factor list
gcnts <- gcnts[ , match(rownames(comparisons), colnames(gcnts))]

#load required packages
pkg1 <- c('pheatmap', 'viridis', 'DESeq2', 'RColorBrewer', 'biomaRt', 'ggrepel')
lapply(pkg1,
       function (x) library(x, character.only = TRUE))

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
load("/Volumes/Gencore/analysis_projects/6368449_Hipschman/deseq_dds_results.rData")


res <- sapply(names(dds),
              function (x) {
                results(dds[[x]], contrast = c(x,2,1))
              },
              simplify = FALSE,
              USE.NAMES = TRUE)

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
         pdf(paste0(x, ".deseq.MDplot.pdf"))
         plotMA(res[[x]], xlim=c(1,1e5), ylim=c(-10,10), main=paste0(x), colSig="grey", colNonSig="grey")
         with(subset(res[[x]], rownames(res[[x]]) %in% rownames(resSig[[x]])), {
          points(resSig[[x]]$baseMean, resSig[[x]]$log2FoldChange, col=ifelse(resSig[[x]]$log2FoldChange>0,"forestgreen","red"), cex = 0.4, lwd = 1, pch=19)
          })
         abline(h = c(-1,1), col = "blue")
         dev.off()
       })

# to set label position: pos = labelPos[[x]],


library(ggplot2)
library(tibble)
library(dplyr)
library(biomaRt)

gcnts <-read.csv("/Volumes/Gencore/sftp/m_fallon/6368449_Hipschman/gene_count.csv", header=TRUE, row.names = "gene_id")

mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))

allGenes <- as.list(getBM(filters = "ensembl_gene_id",
                          attributes = c("ensembl_gene_id", "external_gene_name", "description"),
                          values = rownames(gcnts),
                          mart = mart))
allGenes <- data.frame(allGenes)

results.ordered <- sapply(names(resOrdered),
                          function (x) {
                            x <- as.data.frame(resOrdered[[x]]) %>%
                              rownames_to_column("ENSEMBL") %>%
                              arrange(padj)
                            },
                          simplify = FALSE,
                          USE.NAMES = TRUE)

res.ordered.names <- lapply(results.ordered,
                            function (x) {
                              if (dim(x)[[1]] > 0) {
                                merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                              }
                            })

lapply(names(res.ordered.names),
       function (x) {
         res.ordered.names[[x]] %>%
           ggplot(aes(x = res.ordered.names[[x]]$log2FoldChange,
                      y = -log10(res.ordered.names[[x]]$padj))) + 
           geom_point() + geom_abline(intercept = 2, slope = 0) +
           xlab("log2FoldChange") + ylab("-log10 padj") + ggtitle(paste0(x))
         ggsave(paste0(x, ".deseq.volcano.png"))
         })


pdf("allgenes.deseq.heatmap.pdf")
pheatmap(norm_counts[[1]],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts",
         show_rownames = FALSE)
dev.off()

#### Fibrosis Genes ####
setwd("/Volumes/Gencore/analysis_projects/6368449_Hipschman")

gois <- unique(read.delim("fibrosis.genelist.csv"))
gois <- as.data.frame(gois)

annotated <- merge(gois, allGenes, by.x="GeneName", by.y="external_gene_name")

manual.lines <- data.frame(c("WNT10B", "ENSG00000169884", "wingless-type MMTV integration site family, member 10B [Source:HGNC Symbol;Acc:12775]"),
                           c("ASK1", "ENSG00000197442", "mitogen-activated protein kinase kinase kinase 5 [Source:HGNC Symbol;Acc:HGNC:6857]"),
                           c("COL3A1", "ENSG00000168542", "collagen type III alpha 1 chain [Source:HGNC Symbol;Acc:HGNC:2201]"),
                           c("CTGF", "ENSG00000118523", "connective tissue growth factor [Source:HGNC Symbol;Acc:2500]"),
                           c("NOX2", "ENSG00000165168", "cytochrome b-245 beta chain [Source:HGNC Symbol;Acc:HGNC:2578]"),
                           c("PDGF", "ENSG00000100311", "platelet derived growth factor subunit B [Source:HGNC Symbol;Acc:HGNC:8800]"),
                           c("TIMP-1", "ENSG00000102265", "TIMP metallopeptidase inhibitor 1 [Source:HGNC Symbol;Acc:HGNC:11820]"),
                           c("TIMP-2", "ENSG00000035862", "TIMP metallopeptidase inhibitor 2 [Source:HGNC Symbol;Acc:HGNC:11821]"),
                           c("WNT3A", "ENSG00000154342", "Wnt family member 3A [Source:HGNC Symbol;Acc:HGNC:15983]"))
manual.lines <- t(manual.lines)
rownames(manual.lines) <- NULL
colnames(manual.lines) <- c("GeneName", "ensembl_gene_id", "description")
manual.lines <- as.data.frame(manual.lines)

annotated <- rbind(annotated, manual.lines)

select.degenes <- lapply(res.ordered.names,
                       function (x) {
                         if (dim(x)[[1]] > 0) {
                           merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                         }
                       })

norm.counts <- sapply(names(norm_counts),
                      function (x) {
                        x <- as.data.frame(norm_counts[[x]]) %>%
                          rownames_to_column("ENSEMBL")
                        },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

named.counts <- lapply(norm.counts,
                        function (x) {
                          if (dim(x)[[1]] > 0) {
                            merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                          }
                        })

select.counts <- lapply(named.counts,
                         function (x) {
                           if (dim(x)[[1]] > 0) {
                             merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                           }
                         })

evc1 <- select.degenes$Exp_vs_Control1[select.degenes$Exp_vs_Control1$padj < 0.05,]
evc2 <- select.degenes$Exp_vs_Control2[select.degenes$Exp_vs_Control2$padj < 0.05,]
evcany <- select.degenes$Exp_vs_Controls[select.degenes$Exp_vs_Controls$padj < 0.05,]
c2vc1 <- na.omit(select.degenes$Control2_vs_Control1[select.degenes$Control2_vs_Control1$padj < 0.05,])

write.table(evc1, file = "evc1.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evc2, file = "evc2.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evcany, file = "evcany.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(c2vc1, file = "c2vc1.txt", sep='\t', row.names = FALSE, quote = FALSE)

select.hm <- select.counts$Exp_vs_Control1[c(1:43),]
rownames(select.hm) <- select.hm$GeneName
select.hm <- select.hm[order(select.hm$GeneName),]
select.hm.cat <- as.data.frame(read.delim("fibrosis.geneCats.csv", header = TRUE, sep = ","))
select.hm.cat <- select.hm.cat[order(select.hm.cat$gene_name),]
rownames(select.hm.cat) <- select.hm.cat$gene_name
select.hm.cat <- select.hm.cat[,-1]
#select.hm.cat <- select.counts$Exp_vs_Control1[c(1:43),]
#select.hm.cat <- select.hm.cat[order(select.hm.cat$GeneName),]

ann_colors = list(
  EvC1 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  EvC2 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  C2vC1 = c("no significant difference" = "black", "up" = "green", "down" = "red")
)

# all genes heatmap
pdf("selectgenes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts for Fibrosis Genes",
         labels_row = select.hm$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat$GeneName,
         annotation_colors = ann_colors)
dev.off()

# sig genes heatmap
sig.select.hm.cat <- select.hm.cat %>%
  filter((EvC1 != "no significant difference") | (EvC2 != "no significant difference") | (C2vC1 != "no significant difference"))

sig.select.hm.cat2 <- select.hm.cat[c(5,7,8,12,15,22,24,28:31,33,35,37:39,42,43,44),] 

sig.select.hm <- select.hm[select.hm$GeneName %in% rownames(sig.select.hm.cat),]

rownames(sig.select.hm) <- sig.select.hm$GeneName
pdf("sig.fibrosis.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(sig.select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for Significant Fibrosis Genes",
         labels_row = sig.select.hm$GeneName,
         show_rownames = TRUE,
         annotation_row = sig.select.hm.cat2,
         annotation_colors = ann_colors)
dev.off()

# group average heatmap
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                      Cont1_avg = rowMeans(cbind(Cont1_A, Cont1_B, Cont1_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                      Cont2_avg = rowMeans(cbind(Cont2_A, Cont2_B, Cont2_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                      Exp_avg = rowMeans(cbind(Exp_A, Exp_B, Exp_C)))
select.hm <- as.data.frame(select.hm)
rownames(select.hm) <- select.hm$GeneName
pdf("avg.fibrosis.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,15:17],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for Fibrosis Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# comparative bar chart for GOIs
select.degenes <- lapply(select.degenes,
                         function (x) {
                           x <- x[order(x[['GeneName']]),]
                         })

gene <- c(select.degenes$Exp_vs_Control1$GeneName, 
          select.degenes$Exp_vs_Control2$GeneName,
          select.degenes$Control2_vs_Control1$GeneName)
log2FoldChange <- c(select.degenes$Exp_vs_Control1$log2FoldChange,
                    select.degenes$Exp_vs_Control2$log2FoldChange,
                    select.degenes$Control2_vs_Control1$log2FoldChange)
comparison <- c(rep("Exp vs Cont1", 43), rep("Exp vs Cont2", 43), rep("Cont2 vs Cont1", 43))
select.hm.cat <- rownames_to_column(select.hm.cat, "GeneName")
significance <- merge(select.hm.cat, select.degenes$Exp_vs_Control1, by.x = "GeneName", by.y = "GeneName")
sigcol <- c(significance$EvC1, significance$EvC2, significance$C2vC1)
sigcol <- gsub("no significant difference", " ", sigcol)
sigcol <- gsub("down", "  *", sigcol)
sigcol <- gsub("up", "  *", sigcol)

deData <- data.frame(gene, comparison, log2FoldChange, sigcol)

bars <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, hjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  coord_flip() +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for Fibrosis Genes",
       caption = "Asterisks signify adjusted p-value < 0.05")

ggsave("goiBarPlot.pdf", plot = bars, height = 12, width = 6)

# landscape bar graph
bars_landscape <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, vjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for Fibrosis Genes",
       caption = "Asterisks signify adjusted p-value < 0.05") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggsave("fibrosis.landscape.goiBarPlot.pdf", plot = bars_landscape, height = 4, width = 8)

write.table(file = "goi_normalized_counts.txt", select.counts$Exp_vs_Control1[,-c(11,12)], sep="\t", quote = FALSE)

#### EMT Genes ####

setwd("/Volumes/Gencore/analysis_projects/6368449_Hipschman/EMT")

gois <- unique(read.delim("emt.genelist.csv"))
gois <- as.data.frame(gois)

annotated <- merge(gois, allGenes, by.x="GeneName", by.y="external_gene_name")

manual.lines <- data.frame(c("RNX2", "ENSG00000124813", "RUNX family transcription factor 2 [Source:HGNC Symbol;Acc:HGNC:10472]"))
manual.lines <- t(manual.lines)
rownames(manual.lines) <- NULL
colnames(manual.lines) <- c("GeneName", "ensembl_gene_id", "description")
manual.lines <- as.data.frame(manual.lines)

annotated <- rbind(annotated, manual.lines)

select.degenes <- lapply(res.ordered.names,
                         function (x) {
                           if (dim(x)[[1]] > 0) {
                             merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                           }
                         })

norm.counts <- sapply(names(norm_counts),
                      function (x) {
                        x <- as.data.frame(norm_counts[[x]]) %>%
                          rownames_to_column("ENSEMBL")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

named.counts <- lapply(norm.counts,
                       function (x) {
                         if (dim(x)[[1]] > 0) {
                           merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                         }
                       })

select.counts <- lapply(named.counts,
                        function (x) {
                          if (dim(x)[[1]] > 0) {
                            merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                          }
                        })

evc1 <- select.degenes$Exp_vs_Control1[select.degenes$Exp_vs_Control1$padj < 0.05,]
evc2 <- select.degenes$Exp_vs_Control2[select.degenes$Exp_vs_Control2$padj < 0.05,]
evcany <- select.degenes$Exp_vs_Controls[select.degenes$Exp_vs_Controls$padj < 0.05,]
c2vc1 <- na.omit(select.degenes$Control2_vs_Control1[select.degenes$Control2_vs_Control1$padj < 0.05,])

write.table(evc1, file = "evc1.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evc2, file = "evc2.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evcany, file = "evcany.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(c2vc1, file = "c2vc1.txt", sep='\t', row.names = FALSE, quote = FALSE)

select.hm <- select.counts$Exp_vs_Control1[c(1:11),]
rownames(select.hm) <- select.hm$GeneName
select.hm.cat <- as.data.frame(read.delim("emt.geneCats.csv", header = TRUE, sep = ","))
select.hm.cat <- select.hm.cat[order(select.hm.cat$gene_name),]
rownames(select.hm.cat) <- select.hm.cat$gene_name
select.hm.cat <- select.hm.cat[,-1]

ann_colors = list(
  EvC1 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  EvC2 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  C2vC1 = c("no significant difference" = "black", "up" = "green", "down" = "red")
)

# all genes heatmap
rownames(select.hm) <- select.hm$GeneName
pdf("selectgenes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts for EMT Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# sig genes heatmap
sig.select.hm.cat <- select.hm.cat %>%
  filter((EvC1 != "no significant difference") | (EvC2 != "no significant difference") | (C2vC1 != "no significant difference"))

sig.select.hm <- select.hm[select.hm$GeneName %in% rownames(sig.select.hm.cat),]

rownames(sig.select.hm) <- sig.select.hm$GeneName
pdf("emt.sig.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(sig.select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for EMT Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = sig.select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# group average heatmap
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont1_avg = rowMeans(cbind(Cont1_A, Cont1_B, Cont1_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont2_avg = rowMeans(cbind(Cont2_A, Cont2_B, Cont2_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Exp_avg = rowMeans(cbind(Exp_A, Exp_B, Exp_C)))
select.hm <- as.data.frame(select.hm)
rownames(select.hm) <- select.hm$GeneName
pdf("emt.avg.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,15:17],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for EMT Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# comparative bar chart for GOIs
select.degenes <- lapply(select.degenes,
                         function (x) {
                           x <- x[order(x[['GeneName']]),]
                         })

gene <- c(select.degenes$Exp_vs_Control1$GeneName, 
          select.degenes$Exp_vs_Control2$GeneName,
          select.degenes$Control2_vs_Control1$GeneName)
log2FoldChange <- c(select.degenes$Exp_vs_Control1$log2FoldChange,
                    select.degenes$Exp_vs_Control2$log2FoldChange,
                    select.degenes$Control2_vs_Control1$log2FoldChange)
comparison <- c(rep("Exp vs Cont1", 11), rep("Exp vs Cont2", 11), rep("Cont2 vs Cont1", 11))
select.hm.cat <- rownames_to_column(select.hm.cat, "GeneName")
significance <- merge(select.hm.cat, select.degenes$Exp_vs_Control1, by.x = "GeneName", by.y = "GeneName")
sigcol <- c(significance$EvC1, significance$EvC2, significance$C2vC1)
sigcol <- gsub("no significant difference", " ", sigcol)
sigcol <- gsub("down", "  *", sigcol)
sigcol <- gsub("up", "  *", sigcol)

deData <- data.frame(gene, comparison, log2FoldChange, sigcol)

bars <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, hjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  coord_flip() +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for EMT Genes",
       caption = "Asterisks signify adjusted p-value < 0.05")

ggsave("goiBarPlot.pdf", plot = bars, height = 12, width = 6)

# landscape bar graph
bars_landscape <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, vjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for EMT Genes",
       caption = "Asterisks signify adjusted p-value < 0.05") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggsave("emt.landscape.goiBarPlot.pdf", plot = bars_landscape, height = 4, width = 8)


write.table(file = "goi_normalized_counts.txt", select.counts$Exp_vs_Control1, sep="\t", quote = FALSE)

#### Focal Adhesion Genes ####

setwd("/Volumes/Gencore/analysis_projects/6368449_Hipschman/Focal_Adhesions")

gois <- unique(read.delim("fadh.genelist.csv"))
gois <- as.data.frame(gois)

annotated <- merge(gois, allGenes, by.x="GeneName", by.y="external_gene_name")

manual.lines <- data.frame(c("ILK", "ENSG00000166333", "integrin linked kinase [Source:HGNC Symbol;Acc:HGNC:6040]"))
manual.lines <- t(manual.lines)
rownames(manual.lines) <- NULL
colnames(manual.lines) <- c("GeneName", "ensembl_gene_id", "description")
manual.lines <- as.data.frame(manual.lines)

annotated <- rbind(annotated, manual.lines)

select.degenes <- lapply(res.ordered.names,
                         function (x) {
                           if (dim(x)[[1]] > 0) {
                             merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                           }
                         })

norm.counts <- sapply(names(norm_counts),
                      function (x) {
                        x <- as.data.frame(norm_counts[[x]]) %>%
                          rownames_to_column("ENSEMBL")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

named.counts <- lapply(norm.counts,
                       function (x) {
                         if (dim(x)[[1]] > 0) {
                           merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                         }
                       })

select.counts <- lapply(named.counts,
                        function (x) {
                          if (dim(x)[[1]] > 0) {
                            merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                          }
                        })

evc1 <- select.degenes$Exp_vs_Control1[select.degenes$Exp_vs_Control1$padj < 0.05,]
evc2 <- select.degenes$Exp_vs_Control2[select.degenes$Exp_vs_Control2$padj < 0.05,]
evcany <- select.degenes$Exp_vs_Controls[select.degenes$Exp_vs_Controls$padj < 0.05,]
c2vc1 <- na.omit(select.degenes$Control2_vs_Control1[select.degenes$Control2_vs_Control1$padj < 0.05,])

write.table(evc1, file = "evc1.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evc2, file = "evc2.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evcany, file = "evcany.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(c2vc1, file = "c2vc1.txt", sep='\t', row.names = FALSE, quote = FALSE)

select.hm <- select.counts$Exp_vs_Control1[c(1:30),]
rownames(select.hm) <- select.hm$GeneName
select.hm.cat <- as.data.frame(read.delim("fadh.geneCats.csv", header = TRUE, sep = ","))
select.hm.cat <- select.hm.cat[order(select.hm.cat$gene_name),]
rownames(select.hm.cat) <- select.hm.cat$gene_name
select.hm.cat <- select.hm.cat[,-1]

ann_colors = list(
  EvC1 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  EvC2 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  C2vC1 = c("no significant difference" = "black", "up" = "green", "down" = "red")
)

# all genes heatmap 
rownames(select.hm) <- select.hm$GeneName
pdf("selectgenes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts for Focal Adhesion Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# sig genes heatmap
sig.select.hm.cat <- select.hm.cat %>%
  filter((EvC1 != "no significant difference") | (EvC2 != "no significant difference") | (C2vC1 != "no significant difference"))

sig.select.hm <- select.hm[select.hm$GeneName %in% rownames(sig.select.hm.cat),]

rownames(sig.select.hm) <- sig.select.hm$GeneName
pdf("fadh.sig.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(sig.select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for Focal Adhesion Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = sig.select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# group average heatmap
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont1_avg = rowMeans(cbind(Cont1_A, Cont1_B, Cont1_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont2_avg = rowMeans(cbind(Cont2_A, Cont2_B, Cont2_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Exp_avg = rowMeans(cbind(Exp_A, Exp_B, Exp_C)))
select.hm <- as.data.frame(select.hm)
rownames(select.hm) <- select.hm$GeneName
pdf("fadh.avg.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,15:17],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for Focal Adhesion Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()


# comparative bar chart for GOIs
select.degenes <- lapply(select.degenes,
                         function (x) {
                           x <- x[order(x[['GeneName']]),]
                         })

gene <- c(select.degenes$Exp_vs_Control1$GeneName, 
          select.degenes$Exp_vs_Control2$GeneName,
          select.degenes$Control2_vs_Control1$GeneName)
log2FoldChange <- c(select.degenes$Exp_vs_Control1$log2FoldChange,
                    select.degenes$Exp_vs_Control2$log2FoldChange,
                    select.degenes$Control2_vs_Control1$log2FoldChange)
select.hm.cat <- rownames_to_column(select.hm.cat, "GeneName")
significance <- merge(select.hm.cat, select.degenes$Exp_vs_Control1, by.x = "GeneName", by.y = "GeneName")
comparison <- c(rep("Exp vs Cont1", 30), rep("Exp vs Cont2", 30), rep("Cont2 vs Cont1", 30))
sigcol <- c(significance$EvC1, significance$EvC2, significance$C2vC1)
sigcol <- gsub("no significant difference", " ", sigcol)
sigcol <- gsub("down", "  *", sigcol)
sigcol <- gsub("up", "  *", sigcol)

deData <- data.frame(gene, comparison, log2FoldChange, sigcol)

bars <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, hjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  coord_flip() +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for Focal Adhesion Genes",
       caption = "Asterisks signify adjusted p-value < 0.05")

ggsave("goiBarPlot.pdf", plot = bars, height = 12, width = 6)

# landscape bar graph
bars_landscape <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, vjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for Focal Adhesion Genes",
       caption = "Asterisks signify adjusted p-value < 0.05") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggsave("fadh.landscape.goiBarPlot.pdf", plot = bars_landscape, height = 4, width = 8)


write.table(file = "goi_normalized_counts.txt", select.counts$Exp_vs_Control1, sep="\t", quote = FALSE)

#### HSC Genes ####

setwd("/Volumes/Gencore/analysis_projects/6368449_Hipschman/HSC")

gois <- unique(read.delim("hsc.genelist.csv"))
gois <- as.data.frame(gois)

annotated <- merge(gois, allGenes, by.x="GeneName", by.y="external_gene_name")

manual.lines <- data.frame(c("TNFR1", "ENSG00000067182", "TNF receptor superfamily member 1A [Source:HGNC Symbol;Acc:HGNC:11916]"))
manual.lines <- t(manual.lines)
rownames(manual.lines) <- NULL
colnames(manual.lines) <- c("GeneName", "ensembl_gene_id", "description")
manual.lines <- as.data.frame(manual.lines)

annotated <- rbind(annotated, manual.lines)

select.degenes <- lapply(res.ordered.names,
                         function (x) {
                           if (dim(x)[[1]] > 0) {
                             merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                           }
                         })

norm.counts <- sapply(names(norm_counts),
                      function (x) {
                        x <- as.data.frame(norm_counts[[x]]) %>%
                          rownames_to_column("ENSEMBL")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

named.counts <- lapply(norm.counts,
                       function (x) {
                         if (dim(x)[[1]] > 0) {
                           merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                         }
                       })

select.counts <- lapply(named.counts,
                        function (x) {
                          if (dim(x)[[1]] > 0) {
                            merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                          }
                        })

evc1 <- select.degenes$Exp_vs_Control1[select.degenes$Exp_vs_Control1$padj < 0.05,]
evc2 <- select.degenes$Exp_vs_Control2[select.degenes$Exp_vs_Control2$padj < 0.05,]
evcany <- select.degenes$Exp_vs_Controls[select.degenes$Exp_vs_Controls$padj < 0.05,]
c2vc1 <- na.omit(select.degenes$Control2_vs_Control1[select.degenes$Control2_vs_Control1$padj < 0.05,])

write.table(evc1, file = "evc1.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evc2, file = "evc2.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evcany, file = "evcany.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(c2vc1, file = "c2vc1.txt", sep='\t', row.names = FALSE, quote = FALSE)

select.hm <- select.counts$Exp_vs_Control1[c(1:18),]
rownames(select.hm) <- select.hm$GeneName
select.hm.cat <- as.data.frame(read.delim("hsc.geneCats.csv", header = TRUE, sep = ","))
select.hm.cat <- select.hm.cat[order(select.hm.cat$gene_name),]
rownames(select.hm.cat) <- select.hm.cat$gene_name
select.hm.cat <- select.hm.cat[,-1]

ann_colors = list(
  EvC1 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  EvC2 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  C2vC1 = c("no significant difference" = "black", "up" = "green", "down" = "red")
)

# all genes heatmap
rownames(select.hm) <- select.hm$GeneName
pdf("selectgenes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts for HSC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# sig genes heatmap
sig.select.hm.cat <- select.hm.cat %>%
  filter((EvC1 != "no significant difference") | (EvC2 != "no significant difference") | (C2vC1 != "no significant difference"))

sig.select.hm <- select.hm[select.hm$GeneName %in% rownames(sig.select.hm.cat),]

rownames(sig.select.hm) <- sig.select.hm$GeneName
pdf("hsc.sig.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(sig.select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for HSC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = sig.select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# group average heatmap
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont1_avg = rowMeans(cbind(Cont1_A, Cont1_B, Cont1_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont2_avg = rowMeans(cbind(Cont2_A, Cont2_B, Cont2_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Exp_avg = rowMeans(cbind(Exp_A, Exp_B, Exp_C)))
select.hm <- as.data.frame(select.hm)
rownames(select.hm) <- select.hm$GeneName
pdf("hsc.avg.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,15:17],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for HSC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()


# comparative bar chart for GOIs
select.degenes <- lapply(select.degenes,
                         function (x) {
                           x <- x[order(x[['GeneName']]),]
                         })

gene <- c(select.degenes$Exp_vs_Control1$GeneName, 
          select.degenes$Exp_vs_Control2$GeneName,
          select.degenes$Control2_vs_Control1$GeneName)
log2FoldChange <- c(select.degenes$Exp_vs_Control1$log2FoldChange,
                    select.degenes$Exp_vs_Control2$log2FoldChange,
                    select.degenes$Control2_vs_Control1$log2FoldChange)
select.hm.cat <- rownames_to_column(select.hm.cat, "GeneName")
significance <- merge(select.hm.cat, select.degenes$Exp_vs_Control1, by.x = "GeneName", by.y = "GeneName")
comparison <- c(rep("Exp vs Cont1", 18), rep("Exp vs Cont2", 18), rep("Cont2 vs Cont1", 18))
sigcol <- c(significance$EvC1, significance$EvC2, significance$C2vC1)
sigcol <- gsub("no significant difference", " ", sigcol)
sigcol <- gsub("down", "  *", sigcol)
sigcol <- gsub("up", "  *", sigcol)

deData <- data.frame(gene, comparison, log2FoldChange, sigcol)

bars <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, hjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  coord_flip() +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for HSC Genes",
       caption = "Asterisks signify adjusted p-value < 0.05")

ggsave("goiBarPlot.pdf", plot = bars, height = 12, width = 6)

# landscape bar graph
bars_landscape <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, vjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for HSC Genes",
       caption = "Asterisks signify adjusted p-value < 0.05") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggsave("hsc.landscape.goiBarPlot.pdf", plot = bars_landscape, height = 4, width = 8)



write.table(file = "goi_normalized_counts.txt", select.counts$Exp_vs_Control1, sep="\t", quote = FALSE)

#### HCC Genes ####

setwd("/Volumes/Gencore/analysis_projects/6368449_Hipschman/")

gois <- unique(read.delim("hcc.genelist.csv"))
gois <- as.data.frame(gois)

annotated <- merge(gois, allGenes, by.x="GeneName", by.y="external_gene_name")

manual.lines <- data.frame(c("C-MYC", "ENSG00000136997", "MYC proto-oncogene, bHLH transcription factor [Source:HGNC Symbol;Acc:HGNC:7553]"),
                           c("IGF-1", "ENSG00000017427", "insulin like growth factor 1 [Source:HGNC Symbol;Acc:HGNC:5464]"),
                           c("IGF-1R", "ENSG00000140443", "insulin like growth factor 1 receptor [Source:HGNC Symbol;Acc:HGNC:5465]"),
                           c("IGF-2", "ENSG00000167244", "insulin like growth factor 2 [Source:HGNC Symbol;Acc:HGNC:5466]"),
                           c("IGF-2R", "ENSG00000197081", "insulin like growth factor 2 receptor [Source:HGNC Symbol;Acc:HGNC:5467]"),
                           c("N-RAS", "ENSG00000213281", "NRAS proto-oncogene, GTPase [Source:HGNC Symbol;Acc:HGNC:7989]"),
                           c("VEGFR1", "ENSG00000102755", "fms related receptor tyrosine kinase 1 [Source:HGNC Symbol;Acc:HGNC:3763]"),
                           c("VEGFR2", "ENSG00000128052", "kinase insert domain receptor [Source:HGNC Symbol;Acc:HGNC:6307]"))
manual.lines <- t(manual.lines)
rownames(manual.lines) <- NULL
colnames(manual.lines) <- c("GeneName", "ensembl_gene_id", "description")
manual.lines <- as.data.frame(manual.lines)

annotated <- rbind(annotated, manual.lines)

select.degenes <- lapply(res.ordered.names,
                         function (x) {
                           if (dim(x)[[1]] > 0) {
                             merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                           }
                         })

norm.counts <- sapply(names(norm_counts),
                      function (x) {
                        x <- as.data.frame(norm_counts[[x]]) %>%
                          rownames_to_column("ENSEMBL")
                      },
                      simplify = FALSE,
                      USE.NAMES = TRUE)

named.counts <- lapply(norm.counts,
                       function (x) {
                         if (dim(x)[[1]] > 0) {
                           merge(x, allGenes, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                         }
                       })

select.counts <- lapply(named.counts,
                        function (x) {
                          if (dim(x)[[1]] > 0) {
                            merge(x, annotated, by.x = 'ENSEMBL', by.y = 'ensembl_gene_id')
                          }
                        })

evc1 <- select.degenes$Exp_vs_Control1[select.degenes$Exp_vs_Control1$padj < 0.05,]
evc2 <- select.degenes$Exp_vs_Control2[select.degenes$Exp_vs_Control2$padj < 0.05,]
evcany <- select.degenes$Exp_vs_Controls[select.degenes$Exp_vs_Controls$padj < 0.05,]
c2vc1 <- na.omit(select.degenes$Control2_vs_Control1[select.degenes$Control2_vs_Control1$padj < 0.05,])

write.table(evc1, file = "evc1.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evc2, file = "evc2.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(evcany, file = "evcany.txt", sep='\t', row.names = FALSE, quote = FALSE)
write.table(c2vc1, file = "c2vc1.txt", sep='\t', row.names = FALSE, quote = FALSE)

select.hm <- select.counts$Exp_vs_Control1[c(1:25),]
rownames(select.hm) <- select.hm$GeneName
select.hm.cat <- as.data.frame(read.delim("hcc.geneCats.csv", header = TRUE, sep = ","))
select.hm.cat <- select.hm.cat[order(select.hm.cat$gene_name),]
rownames(select.hm.cat) <- select.hm.cat$gene_name
select.hm.cat <- select.hm.cat[,-1]

ann_colors = list(
  EvC1 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  EvC2 = c("no significant difference" = "black", "up" = "green", "down" = "red"),
  C2vC1 = c("no significant difference" = "black", "up" = "green", "down" = "red")
)

# all genes heatmap

rownames(select.hm) <- select.hm$GeneName
pdf("selectgenes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 6,
         main="Normalized Gene Counts for HCC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# sig genes heatmap
sig.select.hm.cat <- select.hm.cat %>%
  filter((EvC1 != "no significant difference") | (EvC2 != "no significant difference") | (C2vC1 != "no significant difference"))

sig.select.hm <- select.hm[select.hm$GeneName %in% rownames(sig.select.hm.cat),]

rownames(sig.select.hm) <- sig.select.hm$GeneName
pdf("hcc.sig.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(sig.select.hm[,2:10],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for HCC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = sig.select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# group average heatmap
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont1_avg = rowMeans(cbind(Cont1_A, Cont1_B, Cont1_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Cont2_avg = rowMeans(cbind(Cont2_A, Cont2_B, Cont2_C)))
select.hm <- mutate(select.hm %>% 
                      rowwise(), 
                    Exp_avg = rowMeans(cbind(Exp_A, Exp_B, Exp_C)))
select.hm <- as.data.frame(select.hm)
rownames(select.hm) <- select.hm$GeneName
pdf("hcc.avg.genes.deseq.heatmap.pdf", width=10, height=8)
pheatmap(select.hm[,15:17],
         scale="row",
         color=inferno(30),
         fontsize_col = 10,
         main="Normalized Gene Counts for HCC Genes",
         labels_row = select.counts[[1]]$GeneName,
         show_rownames = TRUE,
         annotation_row = select.hm.cat,
         annotation_colors = ann_colors)
dev.off()

# comparative bar chart for GOIs
select.degenes <- lapply(select.degenes,
                         function (x) {
                           x <- x[order(x[['GeneName']]),]
                         })

gene <- c(select.degenes$Exp_vs_Control1$GeneName, 
          select.degenes$Exp_vs_Control2$GeneName,
          select.degenes$Control2_vs_Control1$GeneName)
log2FoldChange <- c(select.degenes$Exp_vs_Control1$log2FoldChange,
                    select.degenes$Exp_vs_Control2$log2FoldChange,
                    select.degenes$Control2_vs_Control1$log2FoldChange)
select.hm.cat <- rownames_to_column(select.hm.cat, "GeneName")
significance <- merge(select.hm.cat, select.degenes$Exp_vs_Control1, by.x = "GeneName", by.y = "GeneName")
comparison <- c(rep("Exp vs Cont1", 25), rep("Exp vs Cont2", 25), rep("Cont2 vs Cont1", 25))
sigcol <- c(significance$EvC1, significance$EvC2, significance$C2vC1)
sigcol <- gsub("no significant difference", " ", sigcol)
sigcol <- gsub("down", "  *", sigcol)
sigcol <- gsub("up", "  *", sigcol)

deData <- data.frame(gene, comparison, log2FoldChange, sigcol)

bars <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, hjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  coord_flip() +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for HCC Genes",
       caption = "Asterisks signify adjusted p-value < 0.05")

ggsave("goiBarPlot.pdf", plot = bars, height = 12, width = 6)

# landscape bar graph
bars_landscape <- ggplot(deData, aes(fill=comparison, y=log2FoldChange, x=gene)) +
  geom_bar(position = position_dodge2(width = 1), stat = "identity") +
  geom_text(aes(label = sigcol, vjust = (0.5 - sign(log2FoldChange)/2)),
            position = position_dodge2(width = 1)) +
  labs(y = "Log 2 Fold Change for Pairwise Comparison",
       x = "Gene",
       title = "Fold Change Difference Between Groups for HCC Genes",
       caption = "Asterisks signify adjusted p-value < 0.05") +
  theme(axis.text.x = element_text(angle = 90, vjust = 1, hjust=1))

ggsave("hcc.landscape.goiBarPlot.pdf", plot = bars_landscape, height = 4, width = 8)


write.table(file = "goi_normalized_counts.txt", select.counts$Exp_vs_Control1, sep="\t", quote = FALSE)

