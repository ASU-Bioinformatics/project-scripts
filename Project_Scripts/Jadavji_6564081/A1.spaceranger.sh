#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -t 2-00:00                    # estimated time needed
#SBATCH --mem=128G
#SBATCH -c 16

#/data/biocore/programs/spaceranger-2.1.1/spaceranger count \
#  --id="Nafisa_MFG-97-40" \
#  --description="Human_Spatial_Tissue_MFG-97-40" \
#  --transcriptome=/data/gencore/databases/spaceranger/refdata-gex-GRCh38-2020-A \
#  --probe-set=/data/gencore/databases/spaceranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
#  --fastqs=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/GZ_Files/MFG-97-40 \
#  --cytaimage=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/Cytassist_picture/assay_CAVG10254_Nafisa_V43A17-371_1694807008_CytAssist/CAVG10254_2023-09-15_19-43-28_Nafisa_V43A17-371_D1_MFG-97-40.tif \
#  --slide=V43A17-371 \
#  --slidefile=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/V43A17-371.gpr \
#  --area=D1 \
#  --reorient-images=true \
#  --localcores=16 \
#  --localmem=128

# MFG-03-41
#--id="Nafisa_MFG-03-41" \
#--description="Human_Spatial_Tissue_MFG-03-41" \
#--fastqs=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/GZ_Files/MFG-03-41 \
#--cytaimage=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/Cytassist_picture/assay_CAVG10254_Nafisa_V43A18-016_1694810271_CytAssist/CAVG10254_2023-09-15_20-37-51_Nafisa_V43A18-016_A1_MFG-03-41.tif \
#--slide=V43A18-016 \
#--slidefile=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/V43A18-016.gpr \
#--area=A1 \

# MFG-10-26
#--id="Nafisa_MFG-10-26" \
#--description="Human_Spatial_Tissue_MFG-10-26" \
#--fastqs=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/GZ_Files/MFG-10-26 \
#--cytaimage=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/Cytassist_picture/assay_CAVG10254_Nafisa_V43A18-016_1694810271_CytAssist/CAVG10254_2023-09-15_20-37-51_Nafisa_V43A18-016_D1_MFG-10-26.tif \
#--slide=V43A18-016 \
#--slidefile=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/V43A18-016.gpr \
#--area=D1 \

# MFG-10-62
/data/biocore/programs/spaceranger-2.1.1/spaceranger count \
  --id="Nafisa_MFG-10-62" \
  --description="Human_Spatial_Tissue_MFG-10-62" \
  --transcriptome=/data/gencore/databases/spaceranger/refdata-gex-GRCh38-2020-A \
  --probe-set=/data/gencore/databases/spaceranger/Visium_Human_Transcriptome_Probe_Set_v2.0_GRCh38-2020-A.csv \
  --fastqs=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/GZ_Files/MFG-10-62 \
  --cytaimage=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/MFG-10-62-cytassist-for-realign.tif \
  --image=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/MFG-10-62_HE-for-realign.TIF \
  --slide=V43A17-371 \
  --slidefile=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/V43A17-371.gpr \
  --area=A1 \
  --loupe-alignment=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/MFG-10-62-for-realign.json \
  --localcores=16 \
  --localmem=128

# MFG-97-40
#--id="Nafisa_MFG-97-40" \
#--description="Human_Spatial_Tissue_MFG-97-40" \
#--fastqs=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/GZ_Files/MFG-97-40 \
#--cytaimage=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/Visium_Analysis/Cytassist_picture/assay_CAVG10254_Nafisa_V43A17-371_1694807008_CytAssist/CAVG10254_2023-09-15_19-43-28_Nafisa_V43A17-371_D1_MFG-97-40.tif \
#--slide=V43A17-371 \
#--slidefile=/data/gencore/analysis_projects/6564081_Jadavji_Spatial/V43A17-371.gpr \
#--area=D1 \
