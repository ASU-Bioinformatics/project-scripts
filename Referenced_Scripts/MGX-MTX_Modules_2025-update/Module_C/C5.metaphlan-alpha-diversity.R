#### LOAD LIBRARIES ####
library(tidyverse)
library(vegan)
library(ggstatsplot)

#### DEFINE VARIABLES AND READ IN DATA ####
metadata <- read.delim("/Volumes/Gencore/analysis_projects/6078853_Otak_RNA/metadata.txt", sep="\t")
metadata <- metadata[order(metadata$Sample.ID),]

setwd("/Volumes/Gencore/analysis_projects/6078853_Otak_RNA/assembly-free-taxonomy/humann-prediction-pipeline/humann3_uniref50/taxonomic-classifications")
bug_df <- read.delim("merged_metaphlan_table_species.txt", sep="\t")

#### CALCULATE PRINCIPAL COMPONENTS ####
colnames(bug_df) <- c("species", metadata$SampleID)

rownames(metadata) <- metadata$SampleID
metadata <- metadata[,-1]

bug_mat = bug_df |> 
  column_to_rownames("species") |> 
  as.matrix() |>
  t()

dist_mat = vegdist(bug_mat)

cmd_res = cmdscale(dist_mat, 
                   k = (nrow(bug_mat) - 1),
                   eig = TRUE)

pcoa_df = tibble(PC1 = cmd_res$points[,1], 
                 PC2 = cmd_res$points[,2])

pcoa_meta = bind_cols(pcoa_df, metadata)

#### PLOT PRINCIPAL COMPONENTS COLORED BY METADATA ####
sc <- viridis::scale_color_viridis("viridis", discrete = FALSE)
plots <- lapply(colnames(metadata),
                function (n) {
                  if (is.character(pcoa_meta[[n]])) {
                    m <- pcoa_meta[[n]]
                    p <- ggplot(pcoa_meta,
                          aes(x=PC1, y=PC2, color=m)) +
                            geom_point(size=3) +
                            labs(color = as.character(n), 
                                  title = paste0("Principal Coordinates by ", n)) +
                            theme_light() +
                            scale_color_brewer(palette = "Set1")
                    ggsave(filename = paste0("pcoa_plot_", n, ".png"), 
                          plot = p,
                          width = 6,
                          height = 4)}
                  else if (is.integer(pcoa_meta[[n]])) {
                    m <- pcoa_meta[[n]]
                    p <- ggplot(pcoa_meta,
                                aes(x=PC1, y=PC2, color=m)) +
                      geom_point(size=3) +
                      sc +
                      labs(color = as.character(n), 
                           title = paste0("Principal Coordinates by ", n)) +
                      theme_light()
                    ggsave(filename = paste0("pcoa_plot_", n, ".png"), 
                           plot = p,
                           width = 6,
                           height = 4)}
                })


wa_data = wascores(cmd_res$points[,1:2], bug_mat) |>
  as_tibble(rownames = 'species')

wa_data

#### ALPHA DIVERSITY ####

categories <- colnames(metadata)
metadata$Shannon <- diversity(bug_mat, index="shannon")
metadata$Simpson <- diversity(bug_mat, index="simpson")
metadata$InverseSimpson <- diversity(bug_mat, index="invsimpson")

write.table(metadata, file="alpha-diversity-with-metadata.txt", quote=FALSE, sep="\t", row.names = FALSE, col.names = TRUE)

metadata$Month <- c("Dec", "Dec", "Feb7", "Feb7", "Feb22", "Feb22")
# the ggbetween stats plots require >=3 samples per category
#metadata <- metadata[,-3]
for (category in categories) {
  #browser()
  if (is.character(metadata[[category]])) {
    print(category)
    plt <- ggbetweenstats(data = metadata,
               x = !!rlang::sym(category),
               y = Simpson,
               pairwise.display = "significant",
               title = paste("Simpson Alpha Diversity for", category))
  
    pdf(paste0("SimpsonAlphaDiversity_", category, ".pdf"))
    print(plt)
    dev.off()
  }}


for (category in categories) {
  if (is.character(metadata[[category]])) {
    plt <- ggbetweenstats(data = metadata,
                        x = !!rlang::sym(category),
                        y = Shannon,
                        pairwise.display = "significant",
                        title = paste("Shannon Alpha Diversity for", category))
    pdf(paste0("ShannonAlphaDiversity_",category,".pdf"))
    print(plt)
    dev.off()
}}



for (category in categories) {
  if (is.character(metadata[[category]])) {
    plt <- ggbetweenstats(data = metadata,
                        x = !!rlang::sym(category),
                        y = InverseSimpson,
                        pairwise.display = "significant",
                        title = paste("Inverse Simpson Alpha Diversity for", category))
    pdf(paste0("InverseSimpsonAlphaDiversity_", category, ".pdf"))
    print(plt)
    dev.off()
}}




#### pathway sandbox ####
setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/assembly-free-taxonomy/humann-prediction-pipeline/humann3_uniref90")

MGX_pwys_raw = read.delim("normalized_pathabundance_stratified.tsv", sep="\t")

setwd("/Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/assembly-free-taxonomy/uniref90-humann-dna")

MTX_pwys_raw = read.delim("normalized_pathabundance_stratified.tsv", sep="\t")

MGX_pwys = MGX_pwys_raw |> 
  separate(X..Pathway, 
           into = c("pwy", "taxa"),
           sep = "\\|") |> 
  filter(!(pwy %in% c("UNMAPPED", "UNINTEGRATED"))) |>
  filter(!(is.na(taxa) | taxa == "unclassified"))

MTX_pwys = MTX_pwys_raw |> 
  separate(X..Pathway, 
           into = c("pwy", "taxa"),
           sep = "\\|") |> 
  filter(!(pwy %in% c("UNMAPPED", "UNINTEGRATED"))) |>
  filter(!(is.na(taxa) | taxa == "unclassified"))

alpha_div = function(pwy_data) {
  
  relative_pwy = pwy_data |>
    mutate(across(c(-pwy, -taxa), ~ .x / sum(.x))) |>
    mutate(across(c(-pwy, -taxa), ~ replace_na(.x, 0))) |>
    select(-taxa) |>
    group_by(pwy)
  
  if (nrow(relative_pwy) == 1) {
    res = relative_pwy |>
      summarise(across(.fns = ~ if_else(sum(.x) > 0,
                                        0,
                                        -1)))
  } else {
    res = relative_pwy |>
      summarise(across(.fns = ~ if_else(sum(.x) > 0,
                                        diversity(.x, index = 'simpson'),
                                        -1)))
  }
  
  sample_values = res |> select(-pwy) |> unlist()
  
  all_null = all(sample_values == -1)
  
  if (all_null) {
    return(NULL)
  } else {
    return(res)
  }
  
}

MGX_a_div = MGX_pwys |> 
  group_split(pwy) |> 
  map(alpha_div) |> 
  bind_rows()

MTX_a_div = MTX_pwys |> 
  group_split(pwy) |> 
  map(alpha_div) |> 
  bind_rows()

threshold = round(0.2 * (ncol(MGX_a_div) - 1))

summary_MGX_a_div = MGX_a_div |> 
  pivot_longer(-pwy, 
               names_to = 'sample_id', 
               values_to = 'alpha_div') |> 
  group_by(pwy) |> 
  mutate(n_detected = sum(alpha_div >= 0)) |> 
  filter(n_detected > threshold) |> 
  summarise(mean_alpha_div = mean(alpha_div[alpha_div >= 0]))

summary_MTX_a_div = MTX_a_div |> 
  pivot_longer(-pwy, 
               names_to = 'sample_id', 
               values_to = 'alpha_div') |> 
  group_by(pwy) |> 
  mutate(n_detected = sum(alpha_div >= 0)) |> 
  filter(n_detected > threshold) |> 
  summarise(mean_alpha_div = mean(alpha_div[alpha_div >= 0]))

mean_alpha_div_df = inner_join(summary_MGX_a_div,
                               summary_MTX_a_div,
                               by = "pwy",
                               suffix = c("_mgx", "_mtx"))


g = ggplot(mean_alpha_div_df, 
           aes(x = mean_alpha_div_mgx, 
               y = mean_alpha_div_mtx)) + 
  geom_point() + 
  theme_bw()

g = g +  
  geom_abline(intercept = 0, slope = 1, color = "grey") + 
  labs(x = "DNA", 
       y = "RNA", 
       title = "Alpha-diversity of Pathway Stratification",
       size = "Log mean\nRNA/DNA ratio") + 
  theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 15), 
        axis.title.y = element_text(size = 15))

g

ggsave(plot = g,
       filename = "outputs/AlphaDiv_DNA_vs_RNA_simple.pdf",
       width = 5,
       height = 5)
