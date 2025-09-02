library("ggplot2")
library("tidyverse")
library("scales")

setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/DM_AIM1/AIM2-panther-overrepresentation")
# the panther input data I'm using is from the merged DEG lists with no log2fc cutoff
# I'm using the Fisher's Exact test with the Bonferroni correction for multiple testing

#### AIM2: c vs b downregulated genes, over-representation of GO Biological Process terms ####
GoTable1 <- read.delim("c_vs_b_down0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$c_vs_b_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$c_vs_b_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$c_vs_b_down0_pantherInput.txt..13.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_c-v-b_down_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: C-Group vs B-Group Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Reduced Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: c vs b upregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("c_vs_b_up0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$c_vs_b_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$c_vs_b_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$c_vs_b_up0_pantherInput.txt..26.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoTable1 <- GoTable1[-3,] #remove unclassified data point
GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment) 

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_c-v-b_up_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, maxCount)) +
  ggtitle("AIM2: c vs b Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Increased Expression in c Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: c vs b upregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("c_vs_b_up0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$c_vs_b_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$c_vs_b_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$c_vs_b_up0_pantherInput.txt..26.
GoTable1 <- GoTable1[-11,] #remove unclassified data point

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_c-v-b_up_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, maxCount)) +
  ggtitle("AIM2: c vs b Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Increased Expression in c Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: c vs b upregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("c_vs_b_up0_goCC-forR.txt", sep="\t")

GoTable1 <- GoTable1[-7,] #remove unclassified data point

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$c_vs_b_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$c_vs_b_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$c_vs_b_up0_pantherInput.txt..26.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_c-v-b_up_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: c vs b Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Increased Expression in c Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()


#### AIM2: calf vs cow downregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("calf_vs_cow_down0_goMF-forR.txt", sep="\t")

GoTable1 <- GoTable1[-7,] #remove unclassified data point

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..10.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_calf-v-cow_down_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, maxCount)) +
  ggtitle("AIM2: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Reduced Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: calf vs cow downregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("calf_vs_cow_down0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..10.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_calf-v-cow_down_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Reduced Expression in Calf Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: calf vs cow downregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("calf_vs_cow_down0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..10.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_calf-v-cow_down_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount)) +
  ggtitle("AIM2: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Reduced Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()


#### AIM2: calf vs cow upregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("calf_vs_cow_up0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_up0_pantherInput.txt..10.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_calf-v-cow_up_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Increased Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un downregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("crypto_vs_un_down0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_down0_pantherInput.txt..676.

GoTable1[GoTable1==" > 100"]<-100
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_down_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Reduced Expression in crypto Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un downregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("crypto_vs_un_down0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_down0_pantherInput.txt..676.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_down_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Reduced Expression in crypto Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un downregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("crypto_vs_un_down0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_down0_pantherInput.txt..676.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_down_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Reduced Expression in crypto Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un upregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("crypto_vs_un_up0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_up0_pantherInput.txt..868.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_up_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Increased Expression in crypto Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un upregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("crypto_vs_un_up0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_up0_pantherInput.txt..868.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_up_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Increased Expression in crypto Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM2: crypto vs un upregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("crypto_vs_un_up0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$crypto_vs_un_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$crypto_vs_un_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$crypto_vs_un_up0_pantherInput.txt..868.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM2_crypto-v-un_up_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0,0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM2: crypto vs un Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Increased Expression in crypto Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()
