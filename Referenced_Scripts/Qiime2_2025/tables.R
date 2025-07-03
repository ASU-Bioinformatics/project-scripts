library(dplyr)
library(knitr)
library(stringr)
library(kableExtra)

setwd("/Volumes/Gencore/analysis_projects/4420250_Marshall_Publication")

#### SHADE ALL DATA ####
statsD <- read.csv("S6-ANCOM-BC-Shade-Stats2digit.csv")

rownames(statsD) <- statsD$Class_ID
statsD <- statsD[, 2:31]

statsS <- lapply(statsD,
                 function(x) {
                   x <- signif(x, 2)
                   x <- ifelse(x==0, x, ifelse(x < 0.01, format(x, scientific = TRUE), format(x, format="fg", scientific = FALSE)))
                   return(x)
                 })

statsL <- lapply(statsS,
                 function(x) {
                   #x <- paste0("$", x)
                   i1 <- grepl("e", x)
                   x[i1] <- paste0(x[i1], "</sup>")
                   x <- sub("e", " x 10<sup>", x)
                   #i2 <- !grepl("10\\^\\{", x)
                   #x[i2] <- paste0(x[i2], "$")
                   return(x)
                 })

statsL <- as.data.frame(statsL)
rownames(statsL) <- rownames(statsD)

statsT <- statsL %>%
  knitr::kable(align = "cccc", format = "html", escape = FALSE,
      col.names = c("SNSP 1", "SNSP 2", "SNSP 3", "SNSP 4",
                    "SUSP 1", "SUSP 2", "SUSP 3", "SUSP 4",
                    "SHPB 1", "SHPB 2", "SHPB 3", "SHPB 4",
                    "Mean SNSP", "Mean SUSP", "Mean SHPB",
                    "Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                    "Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                    "Log2FC", "Standard Error", "p-value", "q-value", "w-score")) %>%
  kable_paper(lightable_options = c("striped", "hover"),
              full_width = F) %>%
  
  column_spec(1:31, width_min = "1in") %>%
  column_spec(1, border_right = T) %>%
  column_spec(16, border_right = T) %>%
  column_spec(21, border_right = T) %>%
  column_spec(26, border_right = T)
  

statsT <- add_header_above(statsT, align = "c",
                           c("Class ID" = 1, "Abundance Values (SNSP = Control, SUSP = Primary, SHPB = Secondary" = 15,
                             "SUSP vs SNSP (Shade Primary vs Shade Control)" = 5, 
                             "SHPB vs SNSP (Shade Secondary vs Shade Control)" = 5,
                             "SHPB vs SUSP (Shade Secondary vs Shade Primary)" = 5))

statsT
save_kable(statsT, "S6-table.html")
save_kable(statsT, "S6-table-format2.html")

#### SHADE ABUNDANCE ONLY ####
statsD <- read.csv("S6-ANCOM-BC-Shade-Abundance-2digit.csv")

rownames(statsD) <- statsD$Class_ID
statsD <- statsD[, 2:16]

statsS <- lapply(statsD,
                 function(x) {
                   x <- signif(x, 2)
                   x <- ifelse(x==0, x, ifelse(x>=0.01, x, format(x, scientific = TRUE)))
                   return(x)
                 })

statsL <- lapply(statsS,
                 function(x) {
                   i1 <- grepl("e", x)
                   x[i1] <- paste0(x[i1], "</sup>")
                   x <- sub("e", " x 10<sup>", x)
                   return(x)
                 })

statsL <- as.data.frame(statsL)
rownames(statsL) <- rownames(statsD)

statsT <- statsL %>%
  knitr::kable(align = "cccc", format = "html", escape = FALSE,
               col.names = c("SNSP 1", "SNSP 2", "SNSP 3", "SNSP 4",
                             "SUSP 1", "SUSP 2", "SUSP 3", "SUSP 4",
                             "SHPB 1", "SHPB 2", "SHPB 3", "SHPB 4",
                             "Mean SNSP", "Mean SUSP", "Mean SHPB")) %>%
  kable_paper(lightable_options = c("striped", "hover"),
              full_width = F) %>%
  column_spec(1:16, width_min = "1in") %>%
  column_spec(1, border_right = T) %>%
  column_spec(5, border_right = T) %>%
  column_spec(9, border_right = T) %>%
  column_spec(13, border_right = T)


statsT <- add_header_above(statsT, align = "c",
                           c("Class ID" = 1,
                             "Shade Control (SNSP)" = 4, 
                             "Shade Primary (SUSP)" = 4,
                             "Shade Secondary (SHPB)" = 4,
                             "Mean Abundance" = 3))

statsT
save_kable(statsT, "shade-abundance-table.html")

#### SHADE STATS ONLY DATA ####
statsD <- read.csv("S6-ANCOM-BC-Shade-StatsOnly-2digit.csv")

rownames(statsD) <- statsD$Class_ID
statsD <- statsD[, 2:16]

statsS <- lapply(statsD,
                 function(x) {
                   x <- signif(x, 2)
                   x <- ifelse(x==0, x, ifelse(x>=0.01, x, format(x, scientific = TRUE)))
                   return(x)
                 })

statsL <- lapply(statsS,
                 function(x) {
                   i1 <- grepl("e", x)
                   x[i1] <- paste0(x[i1], "</sup>")
                   x <- sub("e", " x 10<sup>", x)
                   return(x)
                 })

statsL <- as.data.frame(statsL)
rownames(statsL) <- rownames(statsD)

statsT <- statsL %>%
  knitr::kable(align = "cccc", format = "html", escape = FALSE,
               col.names = c("Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                             "Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                             "Log2FC", "Standard Error", "p-value", "q-value", "w-score")) %>%
  kable_paper(lightable_options = c("striped", "hover"),
              full_width = F) %>%
  
  column_spec(1:16, width_min = "1.2in") %>%
  column_spec(1, border_right = T) %>%
  column_spec(6, border_right = T) %>%
  column_spec(11, border_right = T)


statsT <- add_header_above(statsT, align = "c",
                           c("Class ID" = 1,
                             "SUSP vs SNSP (Shade Primary vs Shade Control)" = 5, 
                             "SHPB vs SNSP (Shade Secondary vs Shade Control)" = 5,
                             "SHPB vs SUSP (Shade Secondary vs Shade Primary)" = 5))

statsT
save_kable(statsT, "shade-stats-table.html")

#### OPEN (SUN) ABUNDANCE ONLY ####
statsD <- read.csv("ANCOM-BC-Sun-Stats.csv")

rownames(statsD) <- statsD$Class_ID
statsD <- statsD[, 2:16]

statsS <- lapply(statsD,
                 function(x) {
                   x <- signif(x, 2)
                   x <- ifelse(x==0, x, ifelse(x>=0.01, x, format(x, scientific = TRUE)))
                   return(x)
                 })

statsL <- lapply(statsS,
                 function(x) {
                   i1 <- grepl("e", x)
                   x[i1] <- paste0(x[i1], "</sup>")
                   x <- sub("e", " x 10<sup>", x)
                   return(x)
                 })

statsL <- as.data.frame(statsL)
rownames(statsL) <- rownames(statsD)

statsT <- statsL %>%
  knitr::kable(align = "cccc", format = "html", escape = FALSE,
               col.names = c("SNSP 1", "SNSP 2", "SNSP 3", "SNSP 4",
                             "SUSP 1", "SUSP 2", "SUSP 3", "SUSP 4",
                             "SHPB 1", "SHPB 2", "SHPB 3", "SHPB 4",
                             "Mean SNSP", "Mean SUSP", "Mean SHPB")) %>%
  kable_paper(lightable_options = c("striped", "hover"),
              full_width = F) %>%
  column_spec(1:16, width_min = "1in") %>%
  column_spec(1, border_right = T) %>%
  column_spec(5, border_right = T) %>%
  column_spec(9, border_right = T) %>%
  column_spec(13, border_right = T)


statsT <- add_header_above(statsT, align = "c",
                           c("Class ID" = 1,
                             "Open Control (SNSPE)" = 4, 
                             "Open Primary (SPE)" = 4,
                             "Open Secondary (SSPB)" = 4,
                             "Mean Abundance" = 3))

statsT
save_kable(statsT, "open-abundance-table.html")

#### OPEN (SUN) STATS ONLY DATA ####
statsD <- read.csv("ANCOM-BC-Sun-Stats.csv")

rownames(statsD) <- statsD$Class_ID
statsD <- statsD[, 2:16]

statsS <- lapply(statsD,
                 function(x) {
                   x <- signif(x, 2)
                   x <- ifelse(x==0, x, ifelse(x>=0.01, x, format(x, scientific = TRUE)))
                   return(x)
                 })

statsL <- lapply(statsS,
                 function(x) {
                   i1 <- grepl("e", x)
                   x[i1] <- paste0(x[i1], "</sup>")
                   x <- sub("e", " x 10<sup>", x)
                   return(x)
                 })

statsL <- as.data.frame(statsL)
rownames(statsL) <- rownames(statsD)

statsT <- statsL %>%
  knitr::kable(align = "cccc", format = "html", escape = FALSE,
               col.names = c("Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                             "Log2FC", "Standard Error", "p-value", "q-value", "w-score",
                             "Log2FC", "Standard Error", "p-value", "q-value", "w-score")) %>%
  kable_paper(lightable_options = c("striped", "hover"),
              full_width = F) %>%
  
  column_spec(1:16, width_min = "1.2in") %>%
  column_spec(1, border_right = T) %>%
  column_spec(6, border_right = T) %>%
  column_spec(11, border_right = T)


statsT <- add_header_above(statsT, align = "c",
                           c("Class ID" = 1,
                             "SPE vs SNPE (Open Primary vs Open Control)" = 5, 
                             "SSPB vs SNPE (Open Secondary vs Open Control)" = 5,
                             "SSPB vs SPE (Open Secondary vs Open Primary)" = 5))

statsT
save_kable(statsT, "open-stats-table.html")
