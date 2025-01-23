if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("gdsfmt")
BiocManager::install("SNPRelate")

library('gdsfmt')
library(SNPRelate)
library(pheatmap)

setwd("/Volumes/Gencore/analysis_projects/7055035_Misra/variant-calls-haploid")

#### SNPRelate ####
gvcf <- "misra.all.genotyped.gatk.g.vcf"
SNPRelate::snpgdsVCF2GDS(gvcf, "misra.all.genotyped.gatk.gds")
SNPRelate::snpgdsSummary("misra.all.genotyped.gatk.gds")
genofile <- SNPRelate::snpgdsOpen("misra.all.genotyped.gatk.gds")
snpset <- SNPRelate::snpgdsLDpruning(genofile, 
                                     ld.threshold=0.8, 
                                     autosome.only = FALSE,
                                     remove.monosnp = FALSE,
                                     maf = NaN,
                                     missing.rate = NaN)
snpset.id <- unlist(unname(snpset))
snpall <- c(1:164)

pca <- SNPRelate::snpgdsPCA(genofile, num.thread=2, autosome.only=FALSE)
pc.percent <- pca$varprop*100
head(round(pc.percent, 2))

subset.pca <- SNPRelate::snpgdsPCA(genofile, snp.i=snpset.id, 
                                   num.thread=2, autosome.only=FALSE)
subset.pc.percent <- subset.pca$varprop*100
head(round(subset.pc.percent, 2))

tab <- data.frame(sample.id = pca$sample.id,
                         EV1 = pca$eigenvect[,1],    # the first eigenvector
                         EV2 = pca$eigenvect[,2],    # the second eigenvector
                         stringsAsFactors = FALSE)
head(tab)

subset.tab <- data.frame(sample.id = subset.pca$sample.id,
                  EV1 = subset.pca$eigenvect[,1],    # the first eigenvector
                  EV2 = subset.pca$eigenvect[,2],    # the second eigenvector
                  stringsAsFactors = FALSE)
head(subset.tab)

lbls <- paste("PC", 1:3, "\n", format(pc.percent[1:3], digits=2), "%", sep="")
tab$col <- c("salmon", "red", "orange", "darkred", "blue", "green", "forestgreen")

pdf("top3-PCA-pairs.pdf", height = 6, width = 6)
pairs(pca$eigenvect[,1:3], col=tab$col, labels=lbls,
      pch=19, main = "Pairwise Visualization of Top Three PCA Dimensions")
dev.off()
pdf("PCA.by.SNPs.pdf", height = 6, width = 6)
plot(tab$EV2, tab$EV1, col=tab$col, pch=19,
     main = "PCA Plot of SNPs per Sample",
     xlab = "PCA Dimension 2", ylab = "PCA Dimension 1")
legend(-0.8, -0.2, legend=c("Mutant 1", "Mutant 2", "Mutant 3", "Mutant 4",
                            "Parent 1", "Parent 2", "Parent 3"),
       col=tab$col, fill=tab$col, cex=0.8)
dev.off()

subset.lbls <- paste("PC", 1:4, "\n", format(pc.percent[1:4], digits=2), "%", sep="")
subset.tab$col <- c("pink", "red", "orange", "yellow", "blue", "green", "forestgreen")
pairs(subset.pca$eigenvect[,1:4], col=subset.tab$col, labels=subset.lbls)
plot(subset.tab$EV2, subset.tab$EV1, xlab="eigenvector 2", ylab="eigenvector 1", col=subset.tab$col, pch=19)

sample.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
ibd <- SNPRelate::snpgdsIBDMoM(genofile, 
                               sample.id=sample.id, 
                               snp.id=snpall,
                               autosome.only = FALSE)

ibd.coeff <- SNPRelate::snpgdsIBDSelection(ibd)
head(ibd.coeff)
ibd.tab <- xtabs(kinship ~ ID1 + ID2, data=ibd.coeff)
ibd.tab[ibd.tab == 0] <- NA
rownames(ibd.tab) <- c("Mutant 1", "Mutant 2", "Mutant 3", "Mutant 4",
                       "Parent 1", "Parent 2")
colnames(ibd.tab) <- c("Mutant 2", "Mutant 3", "Mutant 4",
                       "Parent 1", "Parent 2", "Parent 3")

pdf("ibd-kinship.pdf", height = 8, width = 8)
corrplot(ibd.tab, method="color", is.corr = FALSE, col.lim = c(0.25, 0.5),
         title = "Identity by Descent (Kinship Coefficients)",
         tl.cex = 0.8, tl.col = "darkred", na.label = 'square',
         na.label.col = "white", mar=c(1,0,2,0), type="upper")
dev.off()


plot(ibd.coeff$k0, ibd.coeff$k1, xlim=c(0,1), ylim=c(0,1),
     xlab="k0", ylab="k1", main="All Samples (MoM)")
lines(c(0,1), c(1,0), col="red", lty=2)


ibd.mle <- SNPRelate::snpgdsIBDMLE(genofile, sample.id=sample.id, snp.id=snpall,
                    maf=NaN, missing.rate=NaN, num.thread=2,
                    autosome.only = FALSE, remove.monosnp = FALSE)

ibd.mle.coeff <- SNPRelate::snpgdsIBDSelection(ibd.mle)

plot(ibd.mle.coeff$k0, ibd.mle.coeff$k1, xlim=c(0,1), ylim=c(0,0.000008),
     xlab="k0", ylab="k1", main="IBD for All Samples with All SNPS (MLE)")
lines(c(0,1), c(0.000009,0), col="red", lty=2)

ibs <- SNPRelate::snpgdsIBS(genofile, num.thread=2,
                            autosome.only = FALSE,
                            remove.monosnp = FALSE)
image(ibs$ibs, col=heat.colors(7),
      xlab = "Sample Order: M1, M2, M3, M4, P1, P2, P3",
      ylab = "Sample Order: M1, M2, M3, M4, P1, P2, P3",
      main = "Indentity by State Similarity")

colnames(ibs$ibs) <- c("Mutant 1", "Mutant 2", "Mutant 3", "Mutant 4",
                       "Parent 1", "Parent 2", "Parent 3")
rownames(ibs$ibs) <- c("Mutant 1", "Mutant 2", "Mutant 3", "Mutant 4",
                       "Parent 1", "Parent 2", "Parent 3")
pdf("ibs-similarity.pdf", height = 8, width = 8)
corrplot(ibs$ibs, method = "color", is.corr = FALSE, col.lim = c(0.77, 1),
         title = "Identity by State",
         tl.cex = 0.8, tl.col = "darkred", na.label = 'square',
         na.label.col = "white", mar=c(1,0,2,0), type="upper")
dev.off()

write.table(ibd.coeff, "ibd-table.txt", sep='\t', quote = FALSE)
write.table(ibs$ibs, "ibs-table.txt", sep='\t', quote=FALSE)
