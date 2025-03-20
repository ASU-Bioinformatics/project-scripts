library(dplyr)
library(reshape)
library(gplots)
library(abdiv)
library(ggplot2)
library(ggcorrplot)
library(corrplot)
library(coin)
library(ggsignif)

setwd("/Volumes/Gencore/sftp/t_sandrin/5920921_Keaton/qiime_results/greengenes-classification")

table <- read.delim("gg-rel-level7-table.tsv")
rownames(table) <- table$OTUID
table <- table[,-1]
table.diagnosis <- table %>%
  select(order(unlist(.[1, ])))

diagnosis.cols <- table.diagnosis[1,]
table.diagnosis <- table.diagnosis[-1, ]
table.diagnosis[] <- lapply(table.diagnosis, function(x) {
  as.numeric(x)
})

#### Dissimilarity Indices ####

table.diagnosis <- table.diagnosis %>%
  mutate(Celiac = ifelse(rowSums(table.diagnosis[,1:2]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(CrohnsEst = ifelse(rowSums(table.diagnosis[,3:10]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(CrohnsNaive = ifelse(rowSums(table.diagnosis[,11:12]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(Healthy = ifelse(rowSums(table.diagnosis[,13:16]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(Suspect = ifelse(table.diagnosis[,17] == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(UlcColEst = ifelse(rowSums(table.diagnosis[,18:22]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(UlcColNaive = ifelse(rowSums(table.diagnosis[,23:29]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(UlcPan = ifelse(table.diagnosis[,30] == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(All.Ulcerative.Colitis = ifelse(rowSums(table.diagnosis[,18:29]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(All.Crohns.Disease = ifelse(rowSums(table.diagnosis[,3:12]) == 0, 0, 1))
table.diagnosis <- table.diagnosis %>%
  mutate(All.Disease = ifelse(rowSums(table.diagnosis[,c(1:12,17:30)]) == 0, 0, 1))

jaccard(table.diagnosis$CrohnsEst, table.diagnosis$CrohnsNaive)
jaccard(table.diagnosis$CrohnsNaive, table.diagnosis$Healthy)
jaccard(table.diagnosis$UlcColEst, table.diagnosis$Healthy)
jaccard(table.diagnosis$AllCD, table.diagnosis$Healthy)
jaccard(table.diagnosis$AllUC, table.diagnosis$Healthy)

# distance matrix by diagnosis state
jaccard_distance_matrix.diagnosis <- dist(t(table.diagnosis[,31:38]), method = "binary")
matrix.df.diagnosis <- as.matrix(jaccard_distance_matrix.diagnosis)
matrix.df.diagnosis <- as.data.frame(matrix.df.diagnosis)
write.csv(matrix.df.diagnosis,
          "jacard.presence.absence.beta.diversity.matrix.csv")
png(filename="jacard.presence.absence.beta.diversity.matrix.png", height=6, width=8, units = "in", res = 720)
corrplot(as.matrix(matrix.df.diagnosis), order = "original", type = "upper", col.lim = c(0,1),
         col = COL2('BrBG', 200), tl.cex = 0.8, tl.col = "black",mar=c(0,0,2,0),
         title = "Jaccard Dissimilarity Matrix based on Species Presence/Absence")
dev.off()

# distance matrix by diagnosis group
jaccard_distance_matrix.diagnosis <- dist(t(table.diagnosis[,c(34,39:40)]), method = "binary")
matrix.df.diagnosis <- as.matrix(jaccard_distance_matrix.diagnosis)
matrix.df.diagnosis <- as.data.frame(matrix.df.diagnosis)
write.csv(matrix.df.diagnosis,
          "jacard.presence.absence.beta.diversity.IBDgroup.matrix.csv")
png(filename="jacard.presence.absence.beta.diversity.IBDgroup.matrix.png", height=6, width=8, units = "in", res = 720)
corrplot(as.matrix(matrix.df.diagnosis), order = "original", type = "upper", col.lim = c(0,1), 
         col = COL2('BrBG', 200), tl.cex = 0.8, tl.col = "black",mar=c(0,0,2,0), 
         addCoef.col = "white", addgrid.col = "grey",
         title = "Jaccard Dissimilarity Matrix based on Species Presence/Absence")
dev.off()

# distance matrix by disease state
jaccard_distance_matrix.diagnosis <- dist(t(table.diagnosis[,c(43,34)]), method = "binary")
matrix.df.diagnosis <- as.matrix(jaccard_distance_matrix.diagnosis)
matrix.df.diagnosis <- as.data.frame(matrix.df.diagnosis)
write.csv(matrix.df.diagnosis,
          "jacard.presence.absence.beta.diversity.disease.matrix.csv")
png(filename="jacard.presence.absence.beta.diversity.disease.matrix.png", height=6, width=8, units = "in", res = 720)
corrplot(as.matrix(matrix.df.diagnosis), order = "original", type = "upper", col.lim = c(0,1), 
         col = COL2('BrBG', 200), tl.cex = 0.8, tl.col = "black",mar=c(0,0,2,0), 
         method = "color", addCoef.col = "white", addgrid.col = "grey",
         title = "Jaccard Dissimilarity Matrix based on Species Presence/Absence")
dev.off()

# distance matrix by sample
jaccard_distance_matrix <- dist(t(table.diagnosis[,1:30]), method = "binary")
matrix.df <- as.matrix(jaccard_distance_matrix)
matrix.df <- as.data.frame(matrix.df)
write.csv(matrix.df,
          "jacard.presence.absence.beta.diversity.sample.matrix.csv")
png(filename="jacard.presence.absence.beta.diversity.sample.matrix.png", height=6, width=8, units = "in", res = 720)
corrplot(as.matrix(matrix.df), order = "original", type = "upper", col.lim = c(0,1),
         col = COL2('BrBG', 200), tl.cex = 0.8, tl.col = "black",mar=c(0,0,2,0),
         title = "Jaccard Dissimilarity Matrix based on Species Presence/Absence")
dev.off()

jaclist.CrohnsEst.v.Healthy <- unlist(as.list(matrix.df[3:10,13:16]))
summary(jaclist.CrohnsEst.v.Healthy)
jaclist.CrohnsEst.v.CrohnsNaive <- unlist(as.list(matrix.df[3:10,11:12]))
summary(jaclist.CrohnsEst.v.CrohnsNaive)
jaclist.CrohnsNaive.v.Healthy <- unlist(as.list(matrix.df[11:12,13:16]))
summary(jaclist.CrohnsNaive.v.Healthy)
jaclist.CrohnsEst.v.UlcCol <- unlist(as.list(matrix.df[3:10,18:22]))
summary(jaclist.CrohnsEst.v.UlcCol)
jaclist.CrohnsNaive.v.UlcCol <- unlist(as.list(matrix.df[11:12,18:22]))
summary(jaclist.CrohnsNaive.v.UlcCol)

jaclist.Crohns.v.UlcCol <- unlist(as.list(matrix.df[3:12,18:29]))
jaclist.Crohns.v.Healthy <- unlist(as.list(matrix.df[3:12,13:16]))
jaclist.Healthy.v.UlcCol <- unlist(as.list(matrix.df[13:16,18:29]))
summary(jaclist.Crohns.v.UlcCol)
summary(jaclist.Crohns.v.Healthy)
summary(jaclist.Healthy.v.UlcCol)

# the wilcox and boxplots show that water and surface strains are 
# significantly more similar to each other than either is to the ground strains
# also, that neither group is significantly more similar to the ground strains

wilcox.test(jaclist.Crohns.v.UlcCol, jaclist.Crohns.v.Healthy, exact = FALSE)
wilcox.test(jaclist.Crohns.v.UlcCol, jaclist.Healthy.v.UlcCol, exact = FALSE)
wilcox.test(jaclist.Crohns.v.Healthy, jaclist.Healthy.v.UlcCol, exact = FALSE)

ggplot() + 
  geom_boxplot(aes(x= "All Crohn's v Healthy", 
                   middle = 0.5741, lower = 0.5136, upper = 0.6176, 
                   ymin = 0.4023, ymax = 0.8452), 
               stat = "identity", width = 0.5) +
  geom_jitter(aes(x="All Crohn's v Healthy", y=jaclist.Crohns.v.Healthy),
              size=0.4, width = 0.2) +
  annotate("text", x = 1.5, y = 0.89, label = "p=0.00057", cex=3) +
  annotate("segment", x = 1, xend = 2, y = 0.88, yend = 0.88,
           colour = "blue", linewidth=0.3) +
  geom_boxplot(aes(x= "All Ulcerative Colitis v Healthy", 
                   middle = 0.6372, lower = 0.5954, upper = 0.6912, 
                   ymin = 0.4691, ymax = 0.8265), 
               stat = "identity", width = 0.5) +
  geom_jitter(aes(x="All Ulcerative Colitis v Healthy", y=jaclist.Crohns.v.UlcCol),
              size=0.4, width = 0.2) +
  geom_boxplot(aes(x= "All Crohn's v All Ulcerative Colitis", 
                   middle = 0.6475, lower = 0.5742, upper = 0.7067, 
                   ymin = 0.3953, ymax = 0.8700), 
               stat = "identity", width = 0.5) +
  geom_jitter(aes(x="All Crohn's v All Ulcerative Colitis", y=jaclist.Healthy.v.UlcCol),
              size=0.4, width = 0.2) +
  annotate("text", x = 2.5, y = 0.87, label = "p=0.0015", cex=3) +
  annotate("segment", x = 2, xend = 3, y = 0.86, yend = 0.86,
           colour = "blue", linewidth=0.3) +
  annotate("text", x = 2, y = 0.91, label = "p=0.64", cex=3) +
  annotate("segment", x = 1, xend = 3, y = 0.9, yend = 0.9,
           colour = "red", linewidth=0.3) +
  ggtitle("Jacard Dissimilarity Betwen Groups", subtitle = "By Presence/Absence of Species-Level Taxa from 16S Sequencing") + 
  theme_light() + ylab('Jacard Dissimilarity')
ggsave("Jacard.boxplots.png", height = 6, width = 8)

jaccard_distance_matrix <- dist(t(ann.table.df[,2:15]), method = "binary")
matrix.df <- as.matrix(jaccard_distance_matrix)
matrix.df <- as.data.frame(matrix.df)
write_tsv(matrix.df,
          "jacard.beta.diversity.matrix.tsv",
          col_names = TRUE)
png(filename="jacard.beta.diversity.matrix.png", height=6, width=8, units = "in", res = 720)
corrplot(matrix.df, order = "original", type = "upper", col.lim = c(0,1),
         col = COL2('BrBG', 200), tl.cex = 0.8, tl.col = "black",mar=c(0,0,2,0),
         title = "Jaccard Dissimilarity Matrix based on Gene Presence/Absence")
dev.off()

ann.table.df.noHypo <- ann.table.df[ann.table.df$Annotation != "hypothetical protein",]
ann.table.df.noHypo.uniq <- ann.table.df.noHypo[, c(2:15,17,32:35)] |>
  group_by(Annotation) |>
  summarize(across(.fns = ~ sum(.x))) |>
  mutate(across(-"Annotation", .fns = ~ifelse(.x == 0, 0, 1)))

ann.table.df.noHypo.uniq$ISS.sum <- rowSums(ann.table.df.noHypo.uniq[,2:10])
ann.table.df.noHypo.uniq$ground.sum <- rowSums(ann.table.df.noHypo.uniq[,11:15])
ann.table.df.noHypo.uniq$water.sum <- rowSums(ann.table.df.noHypo.uniq[,4:10])
ann.table.df.noHypo.uniq$surface.sum <- rowSums(ann.table.df.noHypo.uniq[,2:3])
#### Heatmap/Tile Graphs ####
table.ISS <- as.data.frame(ann.table.df[rowSums(ann.table.df[,11:15]) == 0, ])
table.ISS.noHypotheticals <- table.ISS[table.ISS$Annotation != "hypothetical protein",]

table.uniq.values <- ann.table.df[, c(2:15,17,32:35)]
table.uniq.annots <- table.uniq.values |>
  group_by(Annotation) |>
  summarize(across(.fns = ~ sum(.x))) |>
  mutate(across(-"Annotation", .fns = ~ifelse(.x == 0, 0, 1)))

table.uniq.annots.iss <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,11:15]) == 0, ])
table.uniq.annots.iss.long <- melt(table.uniq.annots.iss)
colnames(table.uniq.annots.iss.long) <- c("Annotation", "Strain", "Value")

table.uniq.annots.ground <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,2:10]) == 0, ])

table.uniq.annots.water <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,c(2:3,11:15)]) == 0, ])

table.uniq.annots.surface <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,c(4:15)]) == 0, ])

table.uniq.annots.notAll <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,2:15]) <= 13, ])
table.uniq.annots.notAll.long <- melt(table.uniq.annots.notAll)
colnames(table.uniq.annots.notAll.long) <- c("Annotation", "Strain", "Value")

heatmap.2(as.matrix(table.uniq.annots.notAll[,2:15]),
          cexCol = 0.5)

table.uniq.annots.notTwo <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,2:15]) <= 12, ])
table.uniq.annots.notTwo.long <- melt(table.uniq.annots.notTwo)
colnames(table.uniq.annots.notTwo.long) <- c("Annotation", "Strain", "Value")

heatmap.2(as.matrix(table.uniq.annots.notTwo[,2:15]),
          cexCol = 0.5)

table.uniq.annots.notThree <- as.data.frame(table.uniq.annots[rowSums(table.uniq.annots[,2:15]) <= 11, ])
table.uniq.annots.notThree.long <- melt(table.uniq.annots.notThree)
colnames(table.uniq.annots.notThree.long) <- c("Annotation", "Strain", "Value")

heatmap.2(as.matrix(table.uniq.annots.notThree[,2:15]),
          cexCol = 0.5)
