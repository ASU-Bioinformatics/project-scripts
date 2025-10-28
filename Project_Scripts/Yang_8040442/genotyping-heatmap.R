setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/")

BiocManager::install("ComplexHeatmap")
library(ComplexHeatmap)
library(circlize)

#### E Coli Genotype Heatmaps ####

cft073.gt <- read.delim("/Users/kawoodbu/Documents/TEMP_FOR_SOL/cft073-genotypes-forR.txt")

cft073.mx <- as.matrix(cft073.gt)
cft073.mx.sm <- cft073.mx[rowSums(cft073.mx == ".") < 3, ]
cft073.mx.nz <- cft073.mx[rowSums(cft073.mx == "0") == 0, ]
cft073.mx.nf <- cft073.mx[rowSums(cft073.mx == ".") == 0, ]
cft073.mx.nz.nf <- cft073.mx.nz[rowSums(cft073.mx.nz == ".") == 0, ]
class(cft073.mx.nz.nf) <- "numeric"
cft073.mx.nz.nf <- cft073.mx.nz.nf[rowSums(cft073.mx.nz.nf == "1") < 13, ]

pal <- colorRampPalette(c("white", "darkgrey", "lightgrey"))(3)
colors = structure(pal, names = c("1", "2", "3"))

cft073.hm.nz.nf <- Heatmap(cft073.mx.nz.nf[ , 2:14],
                           col = colors,
                           rect_gp = gpar(col = "darkgrey", lwd = 0.2),
                           column_title = "Genotype Calls for E. Coli SNPs with Confident Calls for All Samples", 
                           column_title_side = "top",
                           column_names_gp = gpar(fontsize = 8),
                           column_labels = c("CAN1102", "X29.F117.AB", "X29.F117.AC",
                                             "X29.F117.AE", "X29.F117.AG", "X29.F117.AI",
                                             "X29.F117.AK", "X29.G117.AB", "X29.G117.AC",
                                             "X29.G117.AD", "X29.G117.AG", 
                                             "X29.G117.AI", "X29.G117.AK"),
                           row_title = "Position in Genome",
                           row_title_gp = gpar(fontsize = 11),
                           row_labels = cft073.mx.nz.nf[ , 1],
                           row_names_gp = gpar(fontsize = 6),
                           name = "genotype",
                           heatmap_legend_param = list(at = 1:3, 
                                                       title = "genotype", 
                                                       legend_gp = gpar(fill = colors), 
                                                       border = "darkgrey"))


png(filename = "cft073.unfiltered.confident-only.snp-genotypes.png",
    width = 8, height = 7, units = "in", res = 720)
draw(cft073.hm.nz.nf)
dev.off()
    
# Filtered
cft073.ft.gt <- read.delim("/Users/kawoodbu/Documents/TEMP_FOR_SOL/cft073-genotypes-filtered-forR.txt")

cft073.ft.mx <- as.matrix(cft073.ft.gt)
cft073.ft.mx.sm <- cft073.ft.mx[rowSums(cft073.ft.mx == ".") < 3, ]
cft073.ft.mx.nz <- cft073.ft.mx[rowSums(cft073.ft.mx == "0") == 0, ]
cft073.ft.mx.nf <- cft073.ft.mx[rowSums(cft073.ft.mx == ".") == 0, ]
cft073.ft.mx.nz.nf <- cft073.ft.mx.nz[rowSums(cft073.ft.mx.nz == ".") == 0, ]
class(cft073.ft.mx.nz.nf) <- "numeric"
cft073.ft.mx.nz.nf <- cft073.ft.mx.nz.nf[rowSums(cft073.ft.mx.nz.nf == "1") < 13, ]

pal <- colorRampPalette(c("white", "darkgrey", "lightgrey"))(3)
colors = structure(pal, names = c("1", "2", "3"))

cft073.ft.hm.nz.nf <- Heatmap(cft073.ft.mx.nz.nf[ , 2:14],
                           col = colors,
                           rect_gp = gpar(col = "lightgrey", lwd = 0.2),
                           column_title = "Genotype Calls for E. Coli SNPs with Hard-Filtered Calls for All Samples", 
                           column_title_side = "top",
                           column_names_gp = gpar(fontsize = 8),
                           column_labels = c("CAN1102", "X29.F117.AB", "X29.F117.AC",
                                             "X29.F117.AE", "X29.F117.AG", "X29.F117.AI",
                                             "X29.F117.AK", "X29.G117.AB", "X29.G117.AC",
                                             "X29.G117.AD", "X29.G117.AG", 
                                             "X29.G117.AI", "X29.G117.AK"),
                           row_title = "Position in Genome",
                           row_title_gp = gpar(fontsize = 11),
                           row_labels = cft073.ft.mx.nz.nf[ , 1],
                           row_names_gp = gpar(fontsize = 7),
                           name = "genotype",
                           heatmap_legend_param = list(at = 1:2, 
                                                       title = "genotype", 
                                                       legend_gp = gpar(fill = colors), 
                                                       border = "darkgrey"))

png(filename = "cft073.filtered.confident-only.snp-genotypes.png",
    width = 8, height = 7, units = "in", res = 720)
draw(cft073.ft.hm.nz.nf)
dev.off()

#### P Aeruginosa Genotype Heatmaps ####

pao1.gt <- read.delim("/Users/kawoodbu/Documents/TEMP_FOR_SOL/pao1-genotypes-unfiltered-forR.txt")

pao1.mx <- as.matrix(pao1.gt)
pao1.mx.sm <- pao1.mx[rowSums(pao1.mx == ".") < 3, ]
pao1.mx.nf <- pao1.mx[rowSums(pao1.mx == ".") == 0, ]
class(pao1.mx.nf) <- "numeric"
pao1.mx.nf <- pao1.mx.nf[rowSums(pao1.mx.nf == "1") < 13, ]

pal <- colorRampPalette(c("white", "darkgrey", "lightgrey"))(3)
colors.pao1 = structure(pal, names = c("0", "1", "."))

pao1.hm.nf <- Heatmap(pao1.mx.nf[ , 2:14],
                           col = colors.pao1,
                           rect_gp = gpar(col = "darkgrey", lwd = 0.2),
                           column_title = "Genotype Calls for P. Aeruginosa SNPs with Confident Calls for All Samples", 
                           column_title_side = "top",
                           column_names_gp = gpar(fontsize = 8),
                           column_labels = c("CAN1103", "X29.F117.AB", "X29.F117.AC",
                                             "X29.F117.AE", "X29.F117.AG", "X29.F117.AI",
                                             "X29.F117.AK", "X29.G117.AB", "X29.G117.AC",
                                             "X29.G117.AD", "X29.G117.AG", 
                                             "X29.G117.AI", "X29.G117.AK"),
                           row_title = "Position in Genome",
                           row_title_gp = gpar(fontsize = 11),
                           row_labels = pao1.mx.nf[ , 1],
                           row_names_gp = gpar(fontsize = 6),
                           name = "genotype",
                           heatmap_legend_param = list(at = 0:3, 
                                                       title = "genotype", 
                                                       legend_gp = gpar(fill = colors.pao1), 
                                                       border = "darkgrey"))


png(filename = "pao1.unfiltered.confident-only.snp-genotypes.png",
    width = 8, height = 7, units = "in", res = 720)
draw(pao1.hm.nf)
dev.off()

# Filtered

# there isn't much point in making an image here, because almost all the SNPs that differ from the reference
# are the same across all the samples in the set.
pao1.ft.gt <- read.delim("/Users/kawoodbu/Documents/TEMP_FOR_SOL/pao1-genotypes-filtered-forR.txt")

pao1.ft.mx <- as.matrix(pao1.ft.gt)
pao1.ft.mx.sm <- pao1.ft.mx[rowSums(pao1.ft.mx == ".") < 3, ]
pao1.ft.mx.nz <- pao1.ft.mx[rowSums(pao1.ft.mx == "0") == 0, ]
pao1.ft.mx.nf <- pao1.ft.mx[rowSums(pao1.ft.mx == ".") == 0, ]
pao1.ft.mx.nz.nf <- pao1.ft.mx.nz[rowSums(pao1.ft.mx.nz == ".") == 0, ]
class(pao1.ft.mx) <- "numeric"
pao1.ft.mx.nz.nf <- pao1.ft.mx.nz.nf[rowSums(pao1.ft.mx.nz.nf == "1") < 13, ]

pao1.ft.hm.nz.nf <- Heatmap(pao1.ft.mx[ , 2:14],
                              col = colors.pao1,
                              rect_gp = gpar(col = "lightgrey", lwd = 0.2),
                              column_title = "Genotype Calls for P. Aeruginosa SNPs with Hard-Filtered Calls for All Samples", 
                              column_title_side = "top",
                              column_names_gp = gpar(fontsize = 8),
                              column_labels = c("CAN1103", "X29.F117.AB", "X29.F117.AC",
                                                "X29.F117.AE", "X29.F117.AG", "X29.F117.AI",
                                                "X29.F117.AK", "X29.G117.AB", "X29.G117.AC",
                                                "X29.G117.AD", "X29.G117.AG", 
                                                "X29.G117.AI", "X29.G117.AK"),
                              row_title = "Position in Genome",
                              row_title_gp = gpar(fontsize = 11),
                              row_labels = pao1.ft.mx[ , 1],
                              row_names_gp = gpar(fontsize = 7),
                              name = "genotype",
                              heatmap_legend_param = list(at = c(0, 1, '.'), 
                                                          title = "genotype", 
                                                          legend_gp = gpar(fill = colors.pao1), 
                                                          border = "darkgrey"))

png(filename = "pao1.filtered.confident-only.snp-genotypes.png",
    width = 8, height = 7, units = "in", res = 720)
draw(pao1.ft.hm.nz.nf)
dev.off()
