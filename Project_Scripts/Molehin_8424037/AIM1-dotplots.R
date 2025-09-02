library("ggplot2")
library("tidyverse")
library("scales")

setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/DM_AIM1/AIM1-panther-overrepresentation")
# the panther input data I'm using is from the merged DEG lists with no log2fc cutoff
# I'm using the Fisher's Exact test with the Bonferroni correction for multiple testing

#### AIM1: calf vs cow downregulated genes, over-representation of GO Molecular Function terms ####
# find a way to create input files from the Panther output instead of typing it manually for the other conditions!
# also to redo this one so the GO IDs are included in the y-axis labels
GoTable1 <- read.delim("calf_vs_cow_down0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..1228.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_down_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Reduced Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: calf vs cow downregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("calf_vs_cow_down0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..1228.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_down_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Reduced Expression in Calf Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: calf vs cow downregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("calf_vs_cow_down0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_down0_pantherInput.txt..1228.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_down_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Reduced Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: calf vs cow upregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("calf_vs_cow_up0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_up0_pantherInput.txt..931.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_up_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Increased Expression in Calf Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: calf vs cow upregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("calf_vs_cow_up0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_up0_pantherInput.txt..931.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_up_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Increased Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: calf vs cow upregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("calf_vs_cow_up0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$calf_vs_cow_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$calf_vs_cow_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$calf_vs_cow_up0_pantherInput.txt..931.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_calf-v-cow_up_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Calf vs Cow Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Increased Expression in Calf Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso downregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("treatments_vs_dmso_down0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..780.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_down_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Reduced Expression in Treatment Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso downregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("treatments_vs_dmso_down0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..780.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_down_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Reduced Expression in Treatment Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso downregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("treatments_vs_dmso_down0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_down0_pantherInput.txt..780.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_down_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Reduced Expression in Treatment Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso upregulated genes, over-representation of GO Molecular Function terms ####
GoTable1 <- read.delim("treatments_vs_dmso_up0_goMF-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..731.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_up_goMF.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Molecular Function Terms in Genes with Increased Expression in Treatment Samples",width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso upregulated genes, over-representation of GO Biological Process terms ####

GoTable1 <- read.delim("treatments_vs_dmso_up0_goBP-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..731.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_up_goBP.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Biological Process Terms in Genes with Increased Expression in Treatment Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### AIM1: treatment vs dmso upregulated genes, over-representation of GO Cellular Component terms ####

GoTable1 <- read.delim("treatments_vs_dmso_up0_goCC-forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.cellular.component.complete
GoTable1$enrichment <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..fold.Enrichment.
GoTable1$pvals <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..P.value.
GoTable1$count <- GoTable1$treatments_vs_dmso_up0_pantherInput.txt..731.
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("AIM1_treatment-v-dmso_up_goCC.png", width = 8, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0.1), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(30)) +
  scale_color_continuous(name = "pvals", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("AIM1: Treatment vs DMSO Functional Over-Representation") +
  labs(subtitle = str_wrap("Top Ten Over-represented GO Cellular Component Terms in Genes with Increased Expression in Treatment Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()
