library("ggplot2")
library("tidyverse")
library("scales")

setwd("/data/gencore/analysis_projects/Jadavji-Spatial-Manuscript/differential-functions-by-status")
# the panther input data I'm using is from the integrated differential features by cluster data from the original Seurat analysis
# I'm using the Fisher's Exact test with the FDR calculation

#### Increased in Control: over-representation of GO BP terms ####
GoTable1 <- read.delim("increased_in_control_goBP_forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$Client.Text.Box.Input..fold.Enrichment.
GoTable1$pvals <- GoTable1$Client.Text.Box.Input..FDR.
GoTable1$count <- GoTable1$Client.Text.Box.Input..228.
GoTable1 <- GoTable1[1:10, ]
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]

# selecting the first 10 before sorting by fold enrichment helps that terms representing clusters of similar terms
# are used, instead of just a lot of terms from a single related cluster

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("increased_in_control_goBP.png", width = 10, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(50)) +
  scale_color_continuous(name = "fdr", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("Functional Terms Overrepresented in Control Patients vs VaD Patients") +
  labs(subtitle = str_wrap("Over-represented GO Biological Process Terms in Genes with Increased Expression in Control Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### Increased in Dementia: over-representation of GO BP terms ####
GoTable1 <- read.delim("increased_in_dementia_goBP_forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.biological.process.complete
GoTable1$enrichment <- GoTable1$Client.Text.Box.Input..fold.Enrichment.
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1$pvals <- GoTable1$Client.Text.Box.Input..FDR.
GoTable1$count <- GoTable1$Client.Text.Box.Input..902.
GoTable1 <- GoTable1[GoTable1$pvals < 0.01, ]
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

# to accomodate the very large number of terms here, I filtered by FDR
# and then sorted by enrichment so I'd grab genes with specificity and impact

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("increased_in_dementia_goBP.png", width = 10, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(50)) +
  scale_color_continuous(name = "fdr", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("Functional Terms Overrepresented in Dementia Patients vs Control Patients") +
  labs(subtitle = str_wrap("Over-represented GO Biological Process Terms in Genes with Increased Expression in Dementia Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### Increased in Dementia: over-representation of GO MF terms ####
GoTable1 <- read.delim("increased_in_dementia_goMF_forR.txt", sep="\t")

GoTable1$GOs <- GoTable1$GO.molecular.function.complete
GoTable1$enrichment <- GoTable1$Client.Text.Box.Input..fold.Enrichment.
GoTable1$enrichment <- as.numeric(GoTable1$enrichment)
GoTable1$pvals <- GoTable1$Client.Text.Box.Input..FDR.
GoTable1$count <- GoTable1$Client.Text.Box.Input..902.
GoTable1 <- GoTable1[GoTable1$pvals < 0.01, ]
GoTable1 <- GoTable1[order(GoTable1$enrichment, decreasing = TRUE),]
GoTable1 <- GoTable1[1:10, ]

# to accomodate the very large number of terms here, I filtered by FDR
# and then sorted by enrichment so I'd grab genes with specificity and impact

GoPlot1 <- ggplot(GoTable1, aes(reorder(GOs, enrichment), enrichment))

minCount = min(GoTable1$count)
meanCount = round(mean(GoTable1$count))
maxCount = max(GoTable1$count)
png("increased_in_dementia_goMF.png", width = 10, height = 6, units = "in", res=720)
GoPlot1 + geom_jitter(position = position_jitter(0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(50)) +
  scale_color_continuous(name = "fdr", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("Functional Terms Overrepresented in Dementia Patients vs Control Patients") +
  labs(subtitle = str_wrap("Over-represented GO Molecular Function Terms in Genes with Increased Expression in Dementia Samples", width = 60),
       y = "Fold Enrichment of GO Term", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

#### Increased in Dementia: over-representation of Panther Pathways ####
PantherTable1 <- read.delim("increased_in_dementia_PantherPathways_forR.txt", sep="\t")

PantherTable1$Panthers <- PantherTable1$PANTHER.Pathways
PantherTable1$enrichment <- PantherTable1$Client.Text.Box.Input..fold.Enrichment.
PantherTable1$enrichment <- as.numeric(PantherTable1$enrichment)
PantherTable1$pvals <- PantherTable1$Client.Text.Box.Input..FDR.
PantherTable1$count <- PantherTable1$Client.Text.Box.Input..902.
PantherTable1 <- PantherTable1[PantherTable1$pvals < 0.01, ]
PantherTable1 <- PantherTable1[order(PantherTable1$enrichment, decreasing = TRUE),]
PantherTable1 <- PantherTable1[1:10, ]

# to accommodate the very large number of terms here, I filtered by FDR
# and then sorted by enrichment so I'd grab genes with specificity and impact

PantherPlot1 <- ggplot(PantherTable1, aes(reorder(Panthers, enrichment), enrichment))

minCount = min(PantherTable1$count)
meanCount = round(mean(PantherTable1$count))
maxCount = max(PantherTable1$count)
png("increased_in_dementia_PantherPathways.png", width = 10, height = 6, units = "in", res=720)
PantherPlot1 + geom_jitter(position = position_jitter(0), aes(color = pvals, size = count)) +
  coord_flip() +
  theme_linedraw() +
  scale_x_discrete(labels = label_wrap(50)) +
  scale_color_continuous(name = "fdr", trans = "log") +
  scale_size_continuous(range = c(4,12),
                        breaks = c(minCount, meanCount, maxCount)) +
  ggtitle("Pathways Overrepresented in Dementia Patients vs Control Patients") +
  labs(subtitle = str_wrap("Over-represented Pathways from PantherDB in Genes with Increased Expression in Dementia Samples", width = 60),
       y = "Fold Enrichment of Pathway", x="") +
  theme(plot.title = element_text(hjust = 0.5),
        plot.subtitle = element_text(hjust = 0.5))
dev.off()

