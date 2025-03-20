library("gplots")
library(tidyverse)
library(vegan)
library(ggstatsplot)
library(reshape)

# the file needed for this script is a tab delimited file with pathways as rownames and metadata values as columns
# the coef value from the Maaslin significant results table is the value for each cell
# only include pathways with a qval > 0.05

setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann50-paired_pipeline/pathways/differential-pathways/MTXmodel-output")

sig.pathways <- read.delim("sig-qval-forR.txt", header = TRUE, row.names = 1, sep = "\t")

path.matrix <- as.matrix(sig.pathways)

melted.paths <- melt(path.matrix)
colnames(melted.paths) <- c("Pathways", "Metadata", "Coefficient")

p <- ggplot(melted.paths, aes(x = Metadata, y = Pathways, fill = Coefficient)) +
  geom_tile() + 
  ggtitle("Differential Pathways by Microorganism", subtitle = "Multivariable Comparative Model") +
  theme_gray(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot")
p

ggsave("maaslin-significant-heatmap.pdf", width = 8, height = 11)  

# exclude groups without biological variation
sig.pathways <- read.delim("sig-qval-forR.txt", header = TRUE, row.names = 1, sep = "\t")
sig.X <- sig.pathways[,!(colnames(sig.pathways) %in% c("Diagnosis.Ulcerative.Pancolitis...Established",
                                                       "Diagnosis.Suspected.IBD..Healthy"))]
sig.X <- sig.X[rowSums(abs(sig.X)) > 0,]

path.matrix.X <- as.matrix(sig.X)

melted.paths.X <- melt(path.matrix.X)
colnames(melted.paths.X) <- c("Pathways", "Metadata", "Coefficient")

p <- ggplot(melted.paths.X, aes(x = Metadata, y = Pathways, fill = Coefficient)) +
  geom_tile() + 
  ggtitle("Differential Pathways by Microorganism", subtitle = "Multivariable Comparative Model") +
  theme_gray(base_size = 5) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
        title = element_text(size = 10),
        plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5),
        plot.title.position = "plot")
p

ggsave("maaslin-significant-heatmap.exclusive.pdf", width = 8, height = 6)  
