# Deconvolution of Spatial RNA Sequencing

## Data Sets:

  1. spatial: data from Nafisa Jadavji and processed by the core

      1. Seurat R objects in an RData file
      2. Spaceranger output data with R scripts available to create Seurat objects

  2. scRNA-seq option 1: data from Darmanis as processed by the Voineageau lab

      1. Darmanis, S. et al. A survey of human brain transcriptome diversity at the single cell level. Proc. Natl Acad. Sci. USA 112, 7285–7290 (2015)
      2. the processed data can be accessed here: https://github.com/VCCRI/CIDR-examples/tree/master/Brain
      3. this dataset is used for the dtangle tutorial, so if CARD doesn't give useful output I'll probably go try dtangle with this dataset following their tutorial.

  3. scRNA-seq option 2 (preferred): data from Lake as processed by the Voineageau lab

      1. Lake, B. B. et al. Integrative single-cell analysis of transcriptional and epigenetic states in the human adult brain. Nat. Biotechnol. 36, 70–80 (2018)
      2. the processed data can be downloaded here: https://github.com/Voineagulab/BrainCellularComposition/blob/main/README.md
      3. the processed data is in a text file in the format of a sparse count matrix, which I know is the needed input for CARD
      4. this dataset provided consistently accurate deconvolution in the study https://www.nature.com/articles/s41467-022-28655-4, so it will be my first choice if the processed data type works.
      5. I cannot find a cell name to cell type annotation table for this dataset so I'm going to try to use data from the Allen Brain Institute, which has this type of data.

  4. scRNA-seq option 3: data from Allen Brain Institute

      1. This data can be downloaded here: https://portal.brain-map.org/atlases-and-data/rnaseq/human-m1-10x
      2.

## Voineageau BrainDeconvShiny App

This is a very user-friendly basic algorithm for quantifying cell type proportions by sample using expression data, and can be accessed here: https://voineagulab.shinyapps.io/BrainDeconvShiny/. It can be used to compare the proportions of cell types between samples, but I don't think it can look at cell type proportions within spots to see the basic spatial pattern across the tissue, or compare gene expression between samples by cell type.

It's still worth running because it is so straightforward, though.

## CARD Algorithm with Lake Reference

A tutorial for the CARD algorithm can be found here: https://yma-lab.github.io/CARD/documentation/04_CARD_Example.html and looks really thorough and detailed. It requires the following input files (the format details are directly taken from the tutorial):

  1. Spatial count data

      1. The spatial transcriptomics count data must be in the format of matrix or sparseMatrix, while each row represents a gene and each column represents a spatial location.
      2. The column names of the spatial data can be in the “XcoordxYcoord” (i.e., 10x10) format, but you can also maintain your original spot names, for example, barcode names.

  2. Spatial location information

      1. The spatial location data must be in the format of data frame while each row represents a spatial location, the first column represents the x coordinate and the second column represents the y coordinate.
      2. The rownames of the spatial location data frame should match exactly with the column names of the spatial_count.
      3. I think I can make this file from the Spaceranger output for each sample, in the file tissue_positions.csv in outs/spatial. That file currently has the barcode, a binary showing whether that spot barcode is present in the image, and the row/column coordinates. I think we just need to pull the barcode columns and the coordinate columns.

  3. scRNAseq (or snRNAseq) count data

      1. The scRNA-seq count data must be in the format of matrix or sparseMatrix, while each row represents a gene and each column represents a cell

  4. scRNAseq (or snRNAseq) metadata

      1. The scRNAseq meta data must be in the format of data frame while each row represents a cell. The rownames of the scRNAseq meta data should match exactly with the column names of the scRNAseq count data.
      2. The sc_meta data must contain the column indicating the cell type assignment for each cell (e.g., “cellType” column in the example sc_meta data).
      3. Sample/subject information should be provided, if there is only one sample, we can add a column by `sc_meta$sampleInfo = "sample1".`
