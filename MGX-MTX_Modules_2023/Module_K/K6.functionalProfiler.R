# the input for this script are the kegg.list files generated from the filtered.annotated.merged.deg files
# (filtered: only DEGs with Kegg terms)
# (annotated: DEG lists combined with annotations)
# (merged: DEGs from all three R programs)

BiocManager::install("MicrobiomeProfiler")

# use the Shiny app to create barplots and dotplots for comparisons with significant Kegg modules
library(MicrobiomeProfiler)
run_MicrobiomeProfiler()

# use this code to generate tables of all results (for any comparison with >2 identified de Kegg terms)
setwd("/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/7.contig-annotation/3.RNA-abundance-profiling-will-need-updating/differential-expression-will-need-updating")

kegglist <- read.table("kegg.list.RNA.filtered_annotated_merged_deg.UlcColNaive.Others.up.txt")

enriched <- enrichKO(kegglist$V1,
                     pvalueCutoff = 0.05,
                     pAdjustMethod = "BH",
                     qvalueCutoff = 0.1)

#setwd("/Volumes/Gencore/sftp/t_sandrin/6209010_Metaomics/7.contig-annotation/3.RNA-abundance-profiling-will-need-updating/differential-expression-will-need-updating")
write_csv(file = "UlcColNaive.Others.up.keggFunctions.table.csv", enriched@result, )
