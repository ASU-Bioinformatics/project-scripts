setwd("/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/5.kraken-assembly-free-analysis/7.bracken-diversity")

library(ggplot2)
library(ggstatsplot)

# copy metadata columns to diversity file prior to uploading
# should be transposed so that SampleID is a column, not a row
alphaDiv <- read.delim("species_alpha-diversity_metadata.txt", header = TRUE, sep = "\t")

metadata <- c("DiseaseState",
              "Diagnosis",
              "Sex",
              "Race")

numericMeta <- c("Age")

for (j in metadata) {
  
  print(j)
  k <- sym(j)
  
  plt.shannon <- ggbetweenstats(
    data = alphaDiv,
    x = !!k,
    y = Shannon,
    pairwise.display = "significant",
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    outlier.label = SampleID)
  
  ggsave(
    width = 7, height = 7,
    paste0("shannon", j, ".png"),
    plt.shannon,
    dpi = 1200
  )
  
  plt.bp <- ggbetweenstats(
    data = alphaDiv,
    x = !!k,
    y = Berger.Parker,
    pairwise.display = "significant",
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    outlier.label = SampleID
  )
  
  ggsave(
    width = 7, height = 7,
    paste0("berger", j, ".png"),
    plt.bp,
    dpi = 1200
  )
  
  plt.invsimp <- ggbetweenstats(
    data = alphaDiv,
    x = !!k,
    y = InverseSimpson,
    pairwise.display = "significant",
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    outlier.label = SampleID
  )
  
  ggsave(
    width = 7, height = 7,
    paste0("inverse.simpson", j, ".png"),
    plt.invsimp,
    dpi = 1200
  )
  
  plt.fisher <- ggbetweenstats(
    data = alphaDiv,
    x = !!k,
    y = Fisher,
    pairwise.display = "significant",
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    outlier.label = SampleID
  )
  
  ggsave(
    width = 7, height = 7,
    paste0("fisher", j, ".png"),
    plt.fisher,
    dpi=1200
  )
   
  plt.simpson <- ggbetweenstats(
    data = alphaDiv,
    x = !!k,
    y = Simpson,
    pairwise.display = "significant",
    type = "nonparametric",
    p.adjust.method = "bonferroni",
    outlier.label = SampleID
  )
  
  ggsave(
    width = 7, height = 7,
    paste0("simpson", j, ".png"),
    plt.simpson,
    dpi=1200
  )
}

for (i in numericMeta) {
  print(i)
  m <- sym(i)
  
  corA <- round(cor(x=alphaDiv[[m]], y=alphaDiv$Shannon), 3)
  shannonPlot <- ggplot(alphaDiv, aes(x=!!m, y=Shannon)) +
    geom_point(size=4,aes(color=Diagnosis)) +
    geom_smooth(method="lm") +
    theme_minimal() + 
    ggtitle(paste0("Shannon Diversity by ", m, ": Correlation = ", corA)) +
    ylab("Shannon Alpha Diversity")
  
  ggsave(
    width = 7, height = 5,
    paste0("shannon-Scatter-", i, ".png"),
    shannonPlot,
    dpi=900,
    bg="white")
  
  corA <- round(cor(x=alphaDiv[[m]], y=alphaDiv$Simpson), 3)
  simpsonPlot <- ggplot(alphaDiv, aes(x=!!m, y=Simpson)) +
    geom_point(size=4,aes(color=Diagnosis)) +
    geom_smooth(method="lm") +
    theme_minimal() + 
    ggtitle(paste0("Simpson Diversity by ", m, ": Correlation = ", corA)) +
    ylab("Simpson Alpha Diversity")
  
  ggsave(
    width = 7, height = 5,
    paste0("simpson-Scatter-", i, ".png"),
    simpsonPlot,
    dpi=900,
    bg="white"
  )
  
  corA <- round(cor(x=alphaDiv[[m]], y=alphaDiv$Berger.Parker), 3)
  bergerPlot <- ggplot(alphaDiv, aes(x=!!m, y=Berger.Parker)) +
    geom_point(size=4,aes(color=Diagnosis)) +
    geom_smooth(method="lm") +
    theme_minimal() + 
    ggtitle(paste0("Berger Diversity by ", m, ": Correlation = ", corA)) +
    ylab("Berger Alpha Diversity")
  
  ggsave(
    width = 7, height = 5,
    paste0("berger-Scatter-", i, ".png"),
    bergerPlot,
    dpi=900,
    bg="white"
  )
  
  corA <- round(cor(x=alphaDiv[[m]], y=alphaDiv$InverseSimpson), 3)
  inverse.simpsonPlot <- ggplot(alphaDiv, aes(x=!!m, y=InverseSimpson)) +
    geom_point(size=4,aes(color=Diagnosis)) +
    geom_smooth(method="lm") +
    theme_minimal() + 
    ggtitle(paste0(" Inverse Simpson Diversity by ", m, ": Correlation = ", corA)) +
    ylab("Inverse Simpson Alpha Diversity")
  
  ggsave(
    width = 7, height = 5,
    paste0("inverse.simpson-Scatter-", i, ".png"),
    inverse.simpsonPlot,
    dpi=900,
    bg="white"
  )
  
  corA <- round(cor(x=alphaDiv[[m]], y=alphaDiv$Fisher), 3)
  fisherPlot <- ggplot(alphaDiv, aes(x=!!m, y=Fisher)) +
    geom_point(size=4,aes(color=Diagnosis)) +
    geom_smooth(method="lm") +
    theme_minimal() + 
    ggtitle(paste0("Fisher Diversity by ", m, ": Correlation = ", corA)) +
    ylab("Fisher Alpha Diversity")
  
  ggsave(
    width = 7, height = 5,
    paste0("fisher-Scatter-", i, ".png"),
    fisherPlot,
    dpi=900,
    bg="white"
  )
}

#### BETA VISUALIZATION ####
library(corrplot)

beta <- read.delim("beta-diversity_forR.txt")

beta <- beta[,-1]
beta = as.data.frame(beta)
beta[beta == "x.xxx"] <- NA
beta <- sapply(beta[,1:30], as.numeric)

colnames(beta) <- c("MB-001", "MB-003", "MB-004", "MB-005", "MB-006", "MB-011",
                    "MB-012", "MB-014", "MB-016", "MB-017", "MB-018", "MB-020",
                    "MB-021", "MB-023", "MB-024", "MB-025", "MB-028", "MB-033",
                    "MB-035", "MB-037", "MB-038", "MB-039", "MB-040", "MB-044",
                    "MB-045", "MB-047", "MB-049", "MB-050", "MB-053", "MB-055")
rownames(beta) <- c("MB-001", "MB-003", "MB-004", "MB-005", "MB-006", "MB-011",
                    "MB-012", "MB-014", "MB-016", "MB-017", "MB-018", "MB-020",
                    "MB-021", "MB-023", "MB-024", "MB-025", "MB-028", "MB-033",
                    "MB-035", "MB-037", "MB-038", "MB-039", "MB-040", "MB-044",
                    "MB-045", "MB-047", "MB-049", "MB-050", "MB-053", "MB-055")

pdf(file = "beta-diversity.pdf", width=6, height=7)
corrplot(beta, method = "color", type = "upper", addgrid.col = "grey",
         is.corr = FALSE, tl.cex=0.8, tl.col = "black")
title(main="Bray-Curtis Dissimilarity Between Samples")
dev.off()
