# prep panther output for R in bash:

# for i in $(find ./ -type f -name "*.txt")
# do
# echo $i
# tail -n +13 $i > "${i%.txt}"-forR.txt
# done

library("ggplot2")
library("tidyverse")
library("scales")

setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod_round2/exp2/rRNA-adjusted")
# the panther input data I'm using is from the DEseq2 output for the model with rRNA abundance adjustment
# including all genes with a p-adjusted value of < 0.05
# I used the Fisher's Exact test with Bonferroni p-value adjustment
# the reference gene list used was the set remaining after filtering and cleaning the dds object

exp1.filelist <- dir(path = ".", recursive = FALSE,
                     pattern = "forR", 
                     full.names = TRUE)

exp1.list <- sapply(exp1.filelist, read.delim, sep="\t", simplify = FALSE, USE.NAMES = TRUE)

exp1.processed <- lapply(exp1.list,
                         function(x) {
                           x$GOs <- x[[1]]
                           x$enrichment <- x[[6]]
                           x$pvals <- x[[7]]
                           x$count <- x[[3]]
                           x[x ==" > 100"]<-100
                           x$enrichment <- as.numeric(x$enrichment)
                           x <- x %>%
                             filter(!grepl('UNCLASSIFIED', GOs))
                           len <- min(15, nrow(x))
                           return(x[order(x$enrichment, decreasing = TRUE), ][1:len, ])
                         })

exp1.plots <- sapply(names(exp1.processed),
                     function(x) {
                       minCount = min(exp1.processed[[x]]$count)
                       meanCount = round(mean(exp1.processed[[x]]$count))
                       maxCount = max(exp1.processed[[x]]$count)
                       plotX <- ggplot(exp1.processed[[x]], aes(reorder(GOs, enrichment), enrichment))
                       return(list("plot" = plotX,
                                   "name" = x,
                                   "min" = minCount,
                                   "mean" = meanCount,
                                   "max" = maxCount,
                                   "length" = length(exp1.processed[[x]]$count)))
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

#split <- str_split((str_split(exp1.plots$`./res.11.v.2.h4.deg.down.sig-slim-goBP-forR.txt`$name, "/")[[1]][2]), "d")
#exp <- substring(split[[1]][1], 1, nchar(split[[1]][1]) - 1)

#split2 <- str_split(split[[1]][2], "\\.")
#ctrl <- substring(split2[[1]][1], 2, nchar(split2[[1]][1]))
#type <- str_split((str_split(exp1.plots$`./res.11.v.2.h4.deg.down.sig-slim-goBP-forR.txt`$name, "/")[[1]][2]), "-")[[1]][3]

#setwd("/Users/kawoodbu/Documents/TEMP_FOR_SOL/8755980_Pallod/Panther_Output/dotplots")
lapply(exp1.plots,
       function(x) {
         split <- str_split((str_split(x$name, "/")[[1]][2]), "d")
         exp <- substring(split[[1]][1], 1, nchar(split[[1]][1]) - 1)
         dir <- str_split((str_split(x$name, "/")[[1]][2]), "\\.")[[1]][7]
         type <- str_split((str_split(x$name, "/")[[1]][2]), "-")[[1]][3]
         if (type == "pathways") {
           name <- "Panther Pathways"
         } else if (type == "goBP") {
           name <- "Biological Process GO Terms"
         } else if (type == "goMF") {
           name <- "Molecular Function GO Terms"
         } else if (type == "goCC") {
           name <- "Cellular Component GO Terms"
         }
         
         if (dir == "up") {
           direction <- "Increased"
         } else {
           direction <- "Reduced"
         }
         
         print(paste0(exp, '.', type, ".png"))
         x$plot + geom_point(aes(reorder(GOs, enrichment), y = enrichment, color = pvals, size = count)) +
           coord_flip() +
           theme_linedraw() +
           scale_x_discrete(labels = label_wrap(30)) +
           scale_color_continuous(name = "pvals", trans = "log") +
           scale_size_continuous(range = c(4,12),
                                 breaks = c(x$min, x$mean, x$max)) +
           ggtitle(paste0("Day 11 Silk with H4_ag vs Day 11 Silk, Functional Over-Representation")) +
           labs(subtitle = str_wrap(paste0("Top Ten Over-represented ", name, " in Genes with ", direction, " Expression in Day 11 Saline"), width = 60),
                y = "Fold Enrichment", x="") +
           theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5))
         height = max(4, 0.75*x$length)
         print(height)
         ggsave(paste0(exp, "_", type, "_", dir, ".png"), width = 8, height = height, units = "in", dpi=720)
       })
