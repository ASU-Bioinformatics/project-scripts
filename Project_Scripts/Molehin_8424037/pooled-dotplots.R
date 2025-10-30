# log in to Phoenix
# interactive -q grp_kawoodbu -p general -t 0-8:00 --mem=32G
# module load r-4.5.1-oj 

# prep panther output for R in bash:

# for i in $(find ./ -type f -name "*.txt")
# do
# echo $i
# tail -n +12 $i > "${i%.txt}"-forR.txt
# done

library("ggplot2")
library("tidyverse")
library("scales")

setwd("/data/gencore/sftp/b_lopez/8424037_Molehin/5a_functional_overrepresentation_pooled/AIM1")
# the panther input data I'm using is from the NOIseq output, including all genes with a probability value >=0.99
# I used the Fisher's Exact test with FDR p-value adjustment

#### AIM1 ####

aim1.filelist <- dir(path = ".", recursive = FALSE,
                            pattern = "forR", 
                            full.names = TRUE)

aim1.list <- sapply(aim1.filelist, read.delim, sep="\t", simplify = FALSE, USE.NAMES = TRUE)

aim1.processed <- lapply(aim1.list,
       function(x) {
         x$GOs <- x[[1]]
         x$enrichment <- x[[6]]
         x$pvals <- x[[8]]
         x$count <- x[[3]]
         len <- min(10, nrow(x))
         x[x ==" > 100"]<-100
         x$enrichment <- as.numeric(x$enrichment) 
         return(x[order(x$enrichment, decreasing = TRUE), ][1:len, ])
       })

aim1.plots <- sapply(names(aim1.processed),
                     function(x) {
                       minCount = min(aim1.processed[[x]]$count)
                       meanCount = round(mean(aim1.processed[[x]]$count))
                       maxCount = max(aim1.processed[[x]]$count)
                       plotX <- ggplot(aim1.processed[[x]], aes(reorder(GOs, enrichment), enrichment))
                       return(list("plot" = plotX,
                                   "name" = x,
                                   "min" = minCount,
                                   "mean" = meanCount,
                                   "max" = maxCount))
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

variable <- "./A_Cow_Pam3csk4_v_A_Cow_DMSO_up_PantherPathways-forR.txt"
exp <- str_split((str_split(variable, "/")[[1]][2]), "_")

lapply(aim1.plots,
       function(x) {
         split <- str_split((str_split(x$name, "/")[[1]][2]), "_")
         exp <- paste0(split[[1]][1], "-", split[[1]][2], "-", split[[1]][3])
         ctrl <- paste0(split[[1]][5], "-", split[[1]][6], "-", split[[1]][7])
         type <- str_split(split[[1]][9], "-")[[1]][1]
         if (type == "PantherPathways") {
           name <- "Panther Pathways"
         } else if (type == "goBP") {
           name <- "Biological Process GO Terms"
         } else if (type == "goMF") {
           name <- "Molecular Function GO Terms"
         } else if (type == "goCC") {
           name <- "Cellular Component GO Terms"
         }
         if (split[[1]][8] == "up") {
           direction <- "Increased"
         } else {
           direction <- "Reduced"
         }
         print(paste0(exp,"_v_", ctrl, "_", type, ".png"))
         #png(paste0(exp,"_v_", ctrl, "_", type, ".png"), width = 8, height = 6, units = "in", res=720)
         x$plot + geom_point(aes(reorder(GOs, enrichment), y = enrichment, color = pvals, size = count)) +
           coord_flip() +
           theme_linedraw() +
           scale_x_discrete(labels = label_wrap(30)) +
           scale_color_continuous(name = "pvals", trans = "log") +
           scale_size_continuous(range = c(4,12),
                                 breaks = c(x$min, x$mean, x$max)) +
           ggtitle(paste0(exp, " vs ", ctrl, " Functional Over-Representation")) +
           labs(subtitle = str_wrap(paste0("Top Ten Over-represented ", name, " in Genes with ", direction, " Expression in ", exp), width = 60),
               y = "Fold Enrichment", x="") +
           theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5))
         ggsave(paste0(exp,"_v_", ctrl, "_", split[[1]][8], "_", type, ".png"), width = 8, height = 6, units = "in", dpi=720)
       })

#### AIM2 ####

setwd("/data/gencore/sftp/b_lopez/8424037_Molehin/5a_functional_overrepresentation_pooled/AIM2")

aim2.filelist <- dir(path = ".", recursive = FALSE,
                     pattern = "forR", 
                     full.names = TRUE)

aim2.list <- sapply(aim2.filelist, read.delim, sep="\t", simplify = FALSE, USE.NAMES = TRUE)

aim2.processed <- lapply(aim2.list,
                         function(x) {
                           x$GOs <- x[[1]]
                           x$enrichment <- x[[6]]
                           x$pvals <- x[[8]]
                           x$count <- x[[3]]
                           len <- min(10, nrow(x))
                           x[x ==" > 100"]<-100
                           x$enrichment <- as.numeric(x$enrichment) 
                           return(x[order(x$enrichment, decreasing = TRUE), ][1:len, ])
                         })

aim2.plots <- sapply(names(aim2.processed),
                     function(x) {
                       minCount = min(aim2.processed[[x]]$count)
                       meanCount = round(mean(aim2.processed[[x]]$count))
                       maxCount = max(aim2.processed[[x]]$count)
                       plotX <- ggplot(aim2.processed[[x]], aes(reorder(GOs, enrichment), enrichment))
                       return(list("plot" = plotX,
                                   "name" = x,
                                   "min" = minCount,
                                   "mean" = meanCount,
                                   "max" = maxCount))
                     },
                     simplify = FALSE,
                     USE.NAMES = TRUE)

lapply(aim2.plots,
       function(x) {
         split <- str_split((str_split(x$name, "/")[[1]][2]), "_")
         exp <- paste0(split[[1]][1], "-", split[[1]][2], "-", split[[1]][3])
         ctrl <- paste0(split[[1]][5], "-", split[[1]][6], "-", split[[1]][7])
         type <- str_split(split[[1]][9], "-")[[1]][1]
         if (type == "PantherPathways") {
           name <- "Panther Pathways"
         } else if (type == "goBP") {
           name <- "Biological Process GO Terms"
         } else if (type == "goMF") {
           name <- "Molecular Function GO Terms"
         } else if (type == "goCC") {
           name <- "Cellular Component GO Terms"
         }
         if (split[[1]][8] == "up") {
           direction <- "Increased"
         } else {
           direction <- "Reduced"
         }
         print(paste0(exp,"_v_", ctrl, "_", type, ".png"))
         #png(paste0(exp,"_v_", ctrl, "_", type, ".png"), width = 8, height = 6, units = "in", res=720)
         x$plot + geom_point(aes(reorder(GOs, enrichment), y = enrichment, color = pvals, size = count)) +
           coord_flip() +
           theme_linedraw() +
           scale_x_discrete(labels = label_wrap(30)) +
           scale_color_continuous(name = "pvals", trans = "log") +
           scale_size_continuous(range = c(4,12),
                                 breaks = c(x$min, x$mean, x$max)) +
           ggtitle(paste0(exp, " vs ", ctrl, " Functional Over-Representation")) +
           labs(subtitle = str_wrap(paste0("Top Ten Over-represented ", name, " in Genes with ", direction, " Expression in ", exp), width = 60),
                y = "Fold Enrichment", x="") +
           theme(plot.title = element_text(hjust = 0.5),
                 plot.subtitle = element_text(hjust = 0.5))
         ggsave(paste0(exp,"_v_", ctrl, "_", split[[1]][8], "_", type, ".png"), width = 8, height = 6, units = "in", dpi=720)
       })
