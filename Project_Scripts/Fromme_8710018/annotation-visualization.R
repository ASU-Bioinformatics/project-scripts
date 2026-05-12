library('circlize')

setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/chloroplast-assembly/params1")
chloroplast.genome <- read.delim("embplant_pt.K105.complete.graph1.1.path_sequence_renamed.bed",
                                 header = FALSE)
chloroplast.colors <- read.delim("embplant_pt.K105.complete.graph1.1.path_sequence_renamed-extra.bed.txt",
                                 header = FALSE)

chloroplast.genome <- merge(chloroplast.genome, chloroplast.colors)
colnames(chloroplast.genome) <- c("chromosome", "start", "end", "name", "score", "strand", "color")
chloroplast.genome <- unique(chloroplast.genome)

# because inverted repeat A crosses over the join of the circular genome,
# I split it into two sections so the bed file can be processed by circlize
chloroplast.genome[24,] <- c("chr1", 121946, 123150, "inverted repeat A", 1, "-", "black")
chloroplast.genome[128,] <- c("chr1", 1, 7, "", 1, "-", "black")

chloro.plus <- chloroplast.genome[chloroplast.genome$strand == "+", ]
chloro.minus <- chloroplast.genome[chloroplast.genome$strand == "-", ]

chloroplast.genome$start <- as.numeric(chloroplast.genome$start)
chloroplast.genome$end <- as.numeric(chloroplast.genome$end)
chloroplast.genome$score <- as.numeric(chloroplast.genome$score)

chloro.plus$start <- as.numeric(chloro.plus$start)
chloro.plus$end <- as.numeric(chloro.plus$end)
chloro.plus$score <- as.numeric(chloro.plus$score)

chloro.minus$start <- as.numeric(chloro.minus$start)
chloro.minus$end <- as.numeric(chloro.minus$end)
chloro.minus$score <- as.numeric(chloro.minus$score)

major.regions <- data.frame(c("chromosome", "start", "end", "name"),
                            c("chr1", 88768, 89978, "IRb"),
                            c("chr1", 121946, 123150, "IRa"),
                            c("chr1", 1, 7, ""),
                            c("chr1", 8, 88767, "LSC"),
                            c("chr1", 89979, 121945, "SSC"))

png("chloroplast.assembly.annotations.png",
    width = 6, height = 6, units = "in", res = 720)
circos.clear()
circos.initializeCircularGenome( chloroplast.genome, chloroplast.genome$name,
                                 123150, plotType = NULL)

circos.genomicLabels(chloro.plus, labels = chloro.plus$name, facing = "clockwise",
                     track.margin=c(0, 0), connection_height = 0.0001,
                     cex = 0.4, niceFacing = TRUE, side = "outside",
                     padding = 0.001, line_col = "white", col = chloro.plus$color)

circos.track(ylim = c(0,1), panel.fun = function(x,y) {
  chr = chloro.plus$name
  xlim = chloro.plus[, 2:3] 
  ylim = c(0,1)
  circos.genomicRect(chloro.plus[, 2:3], chloro.plus$name, col = chloro.plus$color,
                     ytop = 0, ybottom = 1, track.margin=c(0,0.1,0,0), border = NA)
}, track.height = 0.1, bg.border = NA)

circos.genomicAxis(h = "bottom", direction = "inside")

circos.track(ylim = c(0,1), panel.fun = function(x,y) {
  chr = chloro.minus$name
  xlim = chloro.minus[, 2:3] 
  ylim = c(0,1)
  circos.genomicRect(chloro.minus[, 2:3], chloro.minus$name, col = chloro.minus$color,
                     ytop = 0, ybottom = 0.8, track.margin=c(0,0, 0,0), border = NA)
}, track.height = 0.1, bg.border = NA)

circos.genomicLabels(chloro.minus, labels = chloro.minus$name, facing = "clockwise",
                     track.margin=c(0, 0), connection_height = 0.0001,
                     cex = 0.4, niceFacing = TRUE, side = "inside",
                     padding = 0.001, line_col = "white", col = chloro.minus$color)

legend("center",
       legend = c("photosystem I", "other photosynthesis genes",
                  "inverted repeats", "ribosomal proteins (SSU)",
                  "ribosomal proteins (LSU)", "ribosomal RNAs",
                  "transfer RNAs", "other genes", "hypothetical reading frames"),
       col = c("green", "darkgreen", "black", "red",
               "orange", "blue", "blue4", "purple", "hotpink"),
       pch = 16, pt.cex = 0.8, cex = 0.5, text.col = "black",
       inset = c(0.05, 0.05),
       bty = "n")

dev.off()