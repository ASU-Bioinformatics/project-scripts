library(tidyr)
library(scales)
library(patchwork)
library(magrittr)

#### Abundance Data ####
setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/gg-rel-table3-2024.10")
abund <- read.csv("gg-rel-level3-table-2024.10.csv", header = TRUE)

rowMeans(abund[2:25])

abund$max <- do.call(pmax, abund[2:25])
abund.filt <- abund[abund$max >= 0.01, ]
abund.filt <- abund.filt[grep(";c__", abund.filt$OTU_ID), ]
abund.filt <- abund.filt |>
  separate_wider_delim(OTU_ID, delim = "c__", names = c("upper", "class"))

# get abundance values and means per group for classes only in ANCOM results

setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/ancombc-replicate-snpe-level3")
lfc.v.snpe <- read.csv("lfc_slice.csv")

ancom.abund <- merge(abund, lfc.v.snpe, by.x = "OTU_ID", by.y = "id", all.y = TRUE, all.x = FALSE)
ancom.abund <- ancom.abund[ , c(1:2, 13, 19:25, 3:12, 14:18)]
colnames(ancom.abund) <- c("OTU_ID", "SUSP-1", "SUSP-2", "SUSP-3", "SUSP-4",
                           "SNSP-1", "SNSP-2", "SNSP-3", "SNSP-4",
                           "SPE-1", "SPE-2", "SPE-3", "SPE-4",
                           "SNPE-1", "SNPE-2", "SNPE-3", "SNPE-4",
                           "SHPB-1", "SHPB-2", "SHPB-3", "SHPB-4",
                           "SSPB-1", "SSPB-2", "SSPB-3", "SSPB-4")
ancom.abund$mean.susp <- rowMeans(ancom.abund[2:5])
ancom.abund$mean.snsp <- rowMeans(ancom.abund[6:9])
ancom.abund$mean.spe <- rowMeans(ancom.abund[10:13])
ancom.abund$mean.snpe <- rowMeans(ancom.abund[14:17])
ancom.abund$mean.shpb <- rowMeans(ancom.abund[18:21])
ancom.abund$mean.sspb <- rowMeans(ancom.abund[24:27])

setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/")
write.csv(ancom.abund, file = "ANCOM-BC-Abundance.csv")


#### SUN SAMPLES (SPE vs SNPE, SSPB vs SNPE, SSPB vs SPE) ####
setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/ancombc-replicate-snpe-level3")

lfc.v.snpe <- read.csv("lfc_slice.csv")
qval.v.snpe <- read.csv("q_val_slice.csv")
se.v.snpe <- read.csv("se_slice.csv")

setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/ancombc-replicate-spe-level3")

lfc.v.spe <- read.csv("lfc_slice.csv")
qval.v.spe <- read.csv("q_val_slice.csv")
se.v.spe <- read.csv("se_slice.csv")

lfc.v.snpe <- lfc.v.snpe[grep(";c__", lfc.v.snpe$id), ]
qval.v.snpe <- qval.v.snpe[grep(";c__", qval.v.snpe$id), ]
se.v.snpe <- se.v.snpe[grep(";c__", se.v.snpe$id), ]

lfc.v.spe <- lfc.v.spe[grep(";c__", lfc.v.spe$id), ]
qval.v.spe <- qval.v.spe[grep(";c__", qval.v.spe$id), ]
se.v.spe <- se.v.spe[grep(";c__", se.v.spe$id), ]

lfc.v.snpe <- lfc.v.snpe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
lfc.v.snpe <- lfc.v.snpe[, c(2,4:8)]
colnames(lfc.v.snpe) <- c("class", "shpb.v.snpe.lfc", "snsp.v.snpe.lfc",
                   "spe.v.snpe.lfc", "sspb.v.snpe.lfc", "susp.v.snpe.lfc")

lfc.v.spe <- lfc.v.spe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
lfc.v.spe <- lfc.v.spe[, c(2,4:8)]
colnames(lfc.v.spe) <- c("class", "shpb.v.spe.lfc", "snpe.v.spe.lfc",
                          "snsp.v.spe.lfc", "sspb.v.spe.lfc", "susp.v.spe.lfc")

qval.v.snpe <- qval.v.snpe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
qval.v.snpe <- qval.v.snpe[, c(2,4:8)]
colnames(qval.v.snpe) <- c("class", "shpb.v.snpe.qval", "snsp.v.snpe.qval",
                   "spe.v.snpe.qval", "sspb.v.snpe.qval", "susp.v.snpe.qval")

qval.v.spe <- qval.v.spe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
qval.v.spe <- qval.v.spe[, c(2,4:8)]
colnames(qval.v.spe) <- c("class", "shpb.v.spe.qval", "snpe.v.spe.qval",
                           "snsp.v.spe.qval", "sspb.v.spe.qval", "susp.v.spe.qval")

se.v.snpe <- se.v.snpe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
se.v.snpe <- se.v.snpe[, c(2,4:8)]
colnames(se.v.snpe) <- c("class", "shpb.v.snpe.se", "snsp.v.snpe.se",
                   "spe.v.snpe.se", "sspb.v.snpe.se", "susp.v.snpe.se")

se.v.spe <- se.v.spe |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
se.v.spe <- se.v.spe[, c(2,4:8)]
colnames(se.v.spe) <- c("class", "shpb.v.spe.se", "snpe.v.spe.se",
                         "snsp.v.spe.se", "sspb.v.spe.se", "susp.v.spe.se")

data.snpe <- merge(lfc.v.snpe, qval.v.snpe, by = "class")
data.snpe <- merge(data.snpe, se.v.snpe, by = "class")
data.snpe <- data.snpe[data.snpe$class != "uncultured", ]

data.snpe$shpb.v.snpe.low <- data.snpe$shpb.v.snpe.lfc - 1.96*data.snpe$shpb.v.snpe.se
data.snpe$shpb.v.snpe.high <- data.snpe$shpb.v.snpe.lfc + 1.96*data.snpe$shpb.v.snpe.se
data.snpe$snsp.v.snpe.low <- data.snpe$snsp.v.snpe.lfc - 1.96*data.snpe$snsp.v.snpe.se
data.snpe$snsp.v.snpe.high <- data.snpe$snsp.v.snpe.lfc + 1.96*data.snpe$snsp.v.snpe.se
data.snpe$spe.v.snpe.low <- data.snpe$spe.v.snpe.lfc - 1.96*data.snpe$spe.v.snpe.se
data.snpe$spe.v.snpe.high <- data.snpe$spe.v.snpe.lfc + 1.96*data.snpe$spe.v.snpe.se
data.snpe$sspb.v.snpe.low <- data.snpe$sspb.v.snpe.lfc - 1.96*data.snpe$sspb.v.snpe.se
data.snpe$sspb.v.snpe.high <- data.snpe$sspb.v.snpe.lfc + 1.96*data.snpe$sspb.v.snpe.se
data.snpe$susp.v.snpe.low <- data.snpe$susp.v.snpe.lfc - 1.96*data.snpe$susp.v.snpe.se
data.snpe$susp.v.snpe.high <- data.snpe$susp.v.snpe.lfc + 1.96*data.snpe$susp.v.snpe.se

data.spe <- merge(lfc.v.spe, qval.v.spe, by = "class")
data.spe <- merge(data.spe, se.v.spe, by = "class")
data.spe <- data.spe[data.spe$class != "uncultured", ]

data.spe$shpb.v.spe.low <- data.spe$shpb.v.spe.lfc - 1.96*data.spe$shpb.v.spe.se
data.spe$shpb.v.spe.high <- data.spe$shpb.v.spe.lfc + 1.96*data.spe$shpb.v.spe.se
data.spe$snpe.v.spe.low <- data.spe$snpe.v.spe.lfc - 1.96*data.spe$snpe.v.spe.se
data.spe$snpe.v.spe.high <- data.spe$snpe.v.spe.lfc + 1.96*data.spe$snpe.v.spe.se
data.spe$snsp.v.spe.low <- data.spe$snsp.v.spe.lfc - 1.96*data.spe$snsp.v.spe.se
data.spe$snsp.v.spe.high <- data.spe$snsp.v.spe.lfc + 1.96*data.spe$snsp.v.spe.se
data.spe$sspb.v.spe.low <- data.spe$sspb.v.spe.lfc - 1.96*data.spe$sspb.v.spe.se
data.spe$sspb.v.spe.high <- data.spe$sspb.v.spe.lfc + 1.96*data.spe$sspb.v.spe.se
data.spe$susp.v.spe.low <- data.spe$susp.v.spe.lfc - 1.96*data.spe$susp.v.spe.se
data.spe$susp.v.spe.high <- data.spe$susp.v.spe.lfc + 1.96*data.spe$susp.v.spe.se

data.sun <- merge(data.snpe, data.spe, by = "class")
data.sun <- data.sun[((abs(data.sun$spe.v.snpe.lfc) > 2 | abs(data.sun$sspb.v.snpe.lfc) > 2) | (abs(data.sun$spe.v.snpe.lfc) > 2)), 
                 c(1, 4:5, 9:10, 14:15, 21:24, 30, 35, 40, 48:49)]

data.sun <- data.sun[!((data.sun$spe.v.snpe.qval > 0.05 & data.sun$sspb.v.snpe.qval > 0.05) & !(data.sun$sspb.v.spe.qval > 0.05)), ]
data.sun <- merge(data.sun, abund.filt, by = "class", all = FALSE)[, 1:16]

write.csv(data.shade, file = "test.csv")

p <- 
  data.sun |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p <- p +
  geom_point(aes(x=spe.v.snpe.lfc, fill = spe.v.snpe.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=spe.v.snpe.low, xmax=spe.v.snpe.high)) #+ 
  #xlim(-8, 8)

p <- p +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Primary vs Control", y="")
p

p2 <- 
  data.sun |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p2 <- p2 +
  geom_point(aes(x=sspb.v.snpe.lfc, fill = sspb.v.snpe.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=sspb.v.snpe.low, xmax=sspb.v.snpe.high)) #+
  #xlim(-8,8)

p2 <- p2 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Secondary vs Control", y="")
p2

p3 <- 
  data.sun |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p3 <- p3 +
  geom_point(aes(x=sspb.v.spe.lfc, fill = sspb.v.spe.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=sspb.v.spe.low, xmax=sspb.v.spe.high)) +
  labs(fill = "q-value\n") #+ 
  #xlim(-8, 8)

p3 <- p3 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Secondary vs Primary", y="")
p3

layout <- c(
  area(t = 0, l = 0, b = 30, r = 6), 
  area(t = 0, l = 7, b = 30, r = 12.5),
  area(t = 0, l = 13, b = 30, r = 18.5)
  )
# final plot arrangement
setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2")
patch <- p + p2 + p3 + plot_layout(design = layout)
png("ancom.sun.forest.plot.png", height = 6, width = 10, units = "in", res = 720)
patch + plot_annotation(title = "ANCOM-BC Differential Abundance: Sun Exposure",
                        subtitle = "Log2 fold change for classes comprising at least 1% of the population in one or more samples",
                        caption = "Selected classes have a log2 fold change of at least 2 and a q-value of 0.05 or less in at least one comparison. Bars show 95% confidence intervals.")
dev.off()

#### SHADE SAMPLES (SUSP vs SNSP, SHPB vs SNSP, SHPB vs SUSP) ####
setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/ancombc-replicate-snsp-level3")

lfc.v.snsp <- read.csv("lfc_slice.csv")
qval.v.snsp <- read.csv("q_val_slice.csv")
se.v.snsp <- read.csv("se_slice.csv")

setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2/ancombc-replicate-susp-level3")

lfc.v.susp <- read.csv("lfc_slice.csv")
qval.v.susp <- read.csv("q_val_slice.csv")
se.v.susp <- read.csv("se_slice.csv")

lfc.v.snsp <- lfc.v.snsp[grep(";c__", lfc.v.snsp$id), ]
qval.v.snsp <- qval.v.snsp[grep(";c__", qval.v.snsp$id), ]
se.v.snsp <- se.v.snsp[grep(";c__", se.v.snsp$id), ]

lfc.v.susp <- lfc.v.susp[grep(";c__", lfc.v.susp$id), ]
qval.v.susp <- qval.v.susp[grep(";c__", qval.v.susp$id), ]
se.v.susp <- se.v.susp[grep(";c__", se.v.susp$id), ]

lfc.v.snsp <- lfc.v.snsp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
lfc.v.snsp <- lfc.v.snsp[, c(2,4:8)]
colnames(lfc.v.snsp) <- c("class", "shpb.v.snsp.lfc", "snpe.v.snsp.lfc",
                          "spe.v.snsp.lfc", "sspb.v.snsp.lfc", "susp.v.snsp.lfc")

lfc.v.susp <- lfc.v.susp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
lfc.v.susp <- lfc.v.susp[, c(2,4:8)]
colnames(lfc.v.susp) <- c("class", "shpb.v.susp.lfc", "snpe.v.susp.lfc",
                         "snsp.v.susp.lfc", "spe.v.susp.lfc", "sspb.v.susp.lfc")

qval.v.snsp <- qval.v.snsp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
qval.v.snsp <- qval.v.snsp[, c(2,4:8)]
colnames(qval.v.snsp) <- c("class", "shpb.v.snsp.qval", "snpe.v.snsp.qval",
                           "spe.v.snsp.qval", "sspb.v.snsp.qval", "susp.v.snsp.qval")

qval.v.susp <- qval.v.susp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
qval.v.susp <- qval.v.susp[, c(2,4:8)]
colnames(qval.v.susp) <- c("class", "shpb.v.susp.qval", "snpe.v.susp.qval",
                          "snsp.v.susp.qval", "spe.v.susp.qval", "sspb.v.susp.qval")

se.v.snsp <- se.v.snsp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
se.v.snsp <- se.v.snsp[, c(2,4:8)]
colnames(se.v.snsp) <- c("class", "shpb.v.snsp.se", "snpe.v.snsp.se",
                         "spe.v.snsp.se", "sspb.v.snsp.se", "susp.v.snsp.se")

se.v.susp <- se.v.susp |>
  separate_wider_delim(id, delim = "c__", names = c("upper", "class"))
se.v.susp <- se.v.susp[, c(2,4:8)]
colnames(se.v.susp) <- c("class", "shpb.v.susp.se", "snpe.v.susp.se",
                        "snsp.v.susp.se", "spe.v.susp.se", "sspb.v.susp.se")

data.snsp <- merge(lfc.v.snsp, qval.v.snsp, by = "class")
data.snsp <- merge(data.snsp, se.v.snsp, by = "class")
data.snsp <- data.snsp[data.snsp$class != "uncultured", ]

data.snsp$shpb.v.snsp.low <- data.snsp$shpb.v.snsp.lfc - 1.96*data.snsp$shpb.v.snsp.se
data.snsp$shpb.v.snsp.high <- data.snsp$shpb.v.snsp.lfc + 1.96*data.snsp$shpb.v.snsp.se
data.snsp$snpe.v.snsp.low <- data.snsp$snpe.v.snsp.lfc - 1.96*data.snsp$snpe.v.snsp.se
data.snsp$snpe.v.snsp.high <- data.snsp$snpe.v.snsp.lfc + 1.96*data.snsp$snpe.v.snsp.se
data.snsp$spe.v.snsp.low <- data.snsp$spe.v.snsp.lfc - 1.96*data.snsp$spe.v.snsp.se
data.snsp$spe.v.snsp.high <- data.snsp$spe.v.snsp.lfc + 1.96*data.snsp$spe.v.snsp.se
data.snsp$sspb.v.snsp.low <- data.snsp$sspb.v.snsp.lfc - 1.96*data.snsp$sspb.v.snsp.se
data.snsp$sspb.v.snsp.high <- data.snsp$sspb.v.snsp.lfc + 1.96*data.snsp$sspb.v.snsp.se
data.snsp$susp.v.snsp.low <- data.snsp$susp.v.snsp.lfc - 1.96*data.snsp$susp.v.snsp.se
data.snsp$susp.v.snsp.high <- data.snsp$susp.v.snsp.lfc + 1.96*data.snsp$susp.v.snsp.se

data.susp <- merge(lfc.v.susp, qval.v.susp, by = "class")
data.susp <- merge(data.susp, se.v.susp, by = "class")
data.susp <- data.susp[data.susp$class != "uncultured", ]

data.susp$shpb.v.susp.low <- data.susp$shpb.v.susp.lfc - 1.96*data.susp$shpb.v.susp.se
data.susp$shpb.v.susp.high <- data.susp$shpb.v.susp.lfc + 1.96*data.susp$shpb.v.susp.se
data.susp$snpe.v.susp.low <- data.susp$snpe.v.susp.lfc - 1.96*data.susp$snpe.v.susp.se
data.susp$snpe.v.susp.high <- data.susp$snpe.v.susp.lfc + 1.96*data.susp$snpe.v.susp.se
data.susp$snsp.v.susp.low <- data.susp$snsp.v.susp.lfc - 1.96*data.susp$snsp.v.susp.se
data.susp$snsp.v.susp.high <- data.susp$snsp.v.susp.lfc + 1.96*data.susp$snsp.v.susp.se
data.susp$spe.v.susp.low <- data.susp$spe.v.susp.lfc - 1.96*data.susp$spe.v.susp.se
data.susp$spe.v.susp.high <- data.susp$spe.v.susp.lfc + 1.96*data.susp$spe.v.susp.se
data.susp$sspb.v.susp.low <- data.susp$sspb.v.susp.lfc - 1.96*data.susp$sspb.v.susp.se
data.susp$sspb.v.susp.high <- data.susp$sspb.v.susp.lfc + 1.96*data.susp$sspb.v.susp.se

data.shade <- merge(data.snsp, data.susp, by = "class")
write.csv(data.shade, file = "test.csv")
data.shade <- data.shade[((abs(data.shade$susp.v.snsp.lfc) > 2 | abs(data.shade$shpb.v.snsp.lfc) > 2) | (abs(data.shade$shpb.v.susp.lfc) > 2)), 
                     c(1:2, 6:7, 11:12, 16:18, 25:27, 32, 37, 42:43)]

data.shade <- data.shade[!((data.shade$susp.v.snsp.qval > 0.05 & data.shade$shpb.v.snsp.qval > 0.05) & !(data.shade$shpb.v.susp.qval > 0.05)), ]
data.shade <- merge(data.shade, abund.filt, by = "class", all = FALSE)[, 1:16]

p <- 
  data.shade |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(legend.position = "none") +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p <- p +
  geom_point(aes(x=susp.v.snsp.lfc, fill = susp.v.snsp.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=susp.v.snsp.low, xmax=susp.v.snsp.high))

p <- p +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Primary vs Control", y="")
p

p2 <- 
  data.shade |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  theme(legend.position = "none") +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p2 <- p2 +
  geom_point(aes(x=shpb.v.snsp.lfc, fill = shpb.v.snsp.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=shpb.v.snsp.low, xmax=shpb.v.snsp.high))

p2 <- p2 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Secondary vs Control", y="")
p2

p3 <- 
  data.shade |>
  ggplot(aes(y = fct_rev(class))) + 
  theme_classic() +
  theme(axis.text.y = element_blank()) +
  theme(axis.line.y = element_blank()) +
  theme(axis.ticks.y = element_blank()) +
  scale_y_discrete()  +
  geom_hline(yintercept = seq(1, 22, 6), lwd=0.4, color = "#ABABAB") +
  geom_hline(yintercept = seq(2, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(3, 22, 6), lwd=0.4, color = "#CECECE") +
  geom_hline(yintercept = seq(4, 22, 6), lwd=0.4, color = "white") +
  geom_hline(yintercept = seq(5, 22, 6), lwd=0.4, color = "#EBEBEB") +
  geom_hline(yintercept = seq(6, 22, 6), lwd=0.4, color = "white")

p3 <- p3 +
  geom_point(aes(x=shpb.v.susp.lfc, fill = shpb.v.susp.qval), shape=23, size=3) +
  scale_fill_gradient(low = "purple4", high = "white", limits = c(0, 0.05), oob=squish) +
  geom_linerange(aes(xmin=shpb.v.susp.low, xmax=shpb.v.susp.high)) +
  labs(fill = "q-value\n")

p3 <- p3 +
  geom_vline(xintercept = 0, linetype="dashed") +
  labs(x="Secondary vs Primary", y="")
p3

layout <- c(
  area(t = 0, l = 0, b = 30, r = 6), 
  area(t = 0, l = 7, b = 30, r = 12.5),
  area(t = 0, l = 13, b = 30, r = 18.5)
)
# final plot arrangement
setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication/qiime2")
patch <- p + p2 + p3 + plot_layout(design = layout)
png("ancom.shade.forest.plot.png", height = 6, width = 10, units = "in", res = 720)
patch + plot_annotation(title = "ANCOM-BC Differential Abundance: Shade Exposure",
                        subtitle = "Log2 fold change for classes comprising at least 1% of the population in one or more samples",
                        caption = "Selected classes have a log2 fold change of at least 2 and a q-value of 0.05 or less in at least one comparison. Bars show 95% confidence intervals.")
dev.off()
