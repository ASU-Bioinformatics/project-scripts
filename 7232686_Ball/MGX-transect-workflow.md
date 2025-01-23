# Transect Metagenomic Workflow

## Module A - Prepare Reads

For this module, I'm starting with A1.trim-and-cut from the branch version of these scripts, to take advantage of cutadapt's superior adapter trimming while preserving as many reads as possible.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_A/A1.trim-and-cut.sh \
       /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq

mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/original
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-only
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-unpaired

mv *cut_SCT* /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-only
mv *cut_SQP* /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered
mv *cut_SUN* /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-unpaired
mv *fastq.gz /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/original

```

After adapter sequence and low quality reads are removed, we can move on to assembly-free analysis with Kraken, assembly-free analysis with HUMAnN, and coassembly. Since these are soil samples, we don't need to align to a reference genome to remove host reads.

## Module B - Assembly-Free Analysis with Kraken

I am also using the branch files for this module, because the original files don't have help messages and I'd like to add those and potentially clean up other aspects of the code as I go.

### Script B1

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B1.kraken.sh \
 --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
 --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-kraken \
 --dbDir /data/gencore/databases/kraken/k2_pluspf
```

### Script B2

Here, I'm choosing not to use the levels option because I'm satisfied with the default option (all taxonomic levels).

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B2.bracken.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-kraken/kraken-reports \
  --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bracken \
  --dbDir /data/gencore/databases/kraken/k2_pluspf
```

### Script B3

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B3.bracken-diversity.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bracken/bracken-raws \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bracken/bracken-diversity
```

### Script B4

I am using the default installations of KronaTools and KrakenTools to create the plots with this script, so they're not specified here in the command line call.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B4.krona-plots.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bracken/bracken-reports \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-krona/krona-tables \
        --htmlDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-krona/krona-htmls
```

### Script B5 (in R)

I'm recording the custom parts of the script here in case I need to rerun anything! This code is specifically set up for species-level alpha diversity visualization, but I also ran the other taxonomic levels by replacing all instances of 'species' with the other taxonomic levels names respectively.

```
setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bracken/bracken-diversity/species")
div <- read.delim("species_alpha-diversity.txt")
meta <- read.delim("../../transect-metadata.txt")
sampleids <- c("E1", "E2", "E3", "E4", "E5", "L1", "L2", "L3", "L4", "L5", "M1", "M2", "M3", "M4", "M5")
metadata <- c("Stage",
              "PlantType",
              "Stage_and_PlanType")
# no numeric metadata set; the int values in PlantType are discrete and not a numeric series or continuous variables
```

## Module C - Assembly-Free Analysis with HUMAnN

Again I'm using the branch files here so I can correct or improve any of the code.

### Script C1

This script is the heart of the HUMAnN and MetaPhLAN pipeline. Because the fastq files haven't yet been concatenated, I'm using the default value for `--concat` which ensures that the forward and reverse reads are concatenated into a single file prior to classification. This script will be run twice - once for the uniref50 clustered database, and once for the uniref90 clustered database.

Because the database setting is global, DO NOT start running the u50 script until the u90 script is complete! If running multiple projects at once, DO NOT run u50 for any of them until all u90 runs are complete!

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
  --dbType "u90" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann90

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann50

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq/concatenated \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann50 \
  --list "E1-cut E5-cut L2-cut L5-cut M2-cut M3-cut" --resume --already-concatenated

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq/concatenated \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann50 \
  --list "L2-cut" --resume --already-concatenated
```

### Other C Scripts (2-5)
After this script completes, I run the code in `C2.humann-manipulations.sh` interactively to standardize file formatting and location within the output directory.

Scripts `C3.maaslin.R`, `C4.metaphlan-heatmap.R`, and `C5.metaphlan-alpha-diversity.R` are run interactively in a local R session, allowing figures to be fine-tuned and giving me a look at what's going on with the data.

Since a new version of MaAsLin is available now, I added a new script, `C3A.maaslin3.R` to Module_C and ran it. It allowed me to examine compositional abundance and prevalence comparisons for both pathways and taxa, which is really nice! A lot of custom variables are used in this script, from the input data to the model fit formula, so I'm including the entire R script here.

Just a note to self that I accidentally ran the individual sample parts of C2 before L2 had completed analysis, so I'll need to wait and run those manipulations just for L2 before making the merged files.

Since the MaAsLin3 results do look different for Uniref50 vs UniRef90 datasets, I wanted to include a quote from the developers of HUMAnN about the pros and cons of each:

"In brief, UniRef90 is a good default choice since it is comprehensive, non-redundant, and more likely to contain isofunctional clusters. UniRef50 clusters can be very broad, so there's a risk that the cluster representative might not reflect the function of the homologous sequence(s) found in your dataset. One situation where UniRef50 might be preferable is when dealing with very poorly characterized microbiomes. In that case, requiring reads to map at 90% identity to UniRef90 might be too stringent, and so mapping at 50% identity to UniRef50 could explain a larger portion of sample reads (albeit at reduced resolution)."

```
#### DEFINE CUSTOM VARIABLES ####

# this is only for nonpaired data. Use MTX model scripts for paired metatranscriptomic data (section J)
# to prepare abundance file for R, remove hash from Pathway column and shorten the sample names

metaDir <- "/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/"
metaFile <- paste0(metaDir,"transect-metadata.txt")
projectDir <- "/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann90/"
outputDir <- paste0(projectDir, "pathways/differential-pathways")
inputDir <- paste0(projectDir, "pathways/pathway-abundance/")
taxaDir <- paste0(projectDir, "taxonomic-classifications/")
normalizedAbundance <- paste0(inputDir, "normalized_pathabundance_unstratified.tsv")
metaphlanTaxonomy <- paste0(taxaDir, "merged_metaphlan_table_species.txt")

#### SET DIRECTORY AND UPLOAD METADATA ####

setwd(outputDir)

metadata <- data.frame(read.delim(metaFile, sep="\t"))

rownames(metadata) <- metadata[,1]
metadata <- metadata[,-1]
metadata$reads <- c(32681050, 28668149, 23782422, 33350609, 33083108,
                    33392171, 203558894, 26473730, 26619467, 24251553,
                    25081992, 39481746, 32557645, 23126532, 41927618)

#### ANALYZE MAASLIN ####

# input and format pathway input
pa.input <- data.frame(read.delim(normalizedAbundance, sep="\t"))

rownames(pa.input) <- pa.input[,1]
pa.input <- pa.input[,-1]

pa.input <- as.data.frame(t(pa.input))

# input and format taxonomy input (removed 's_' prefix from taxon names and shortened sample names in text editor)
tax.input <- data.frame(read.delim(metaphlanTaxonomy, sep="\t"))

rownames(tax.input) <- tax.input[,1]
tax.input <- tax.input[,-1]

tax.input <- as.data.frame(t(tax.input))


# factor the categorical variables in metadata
metadata$Stage <-
  factor(metadata$Stage, levels = c('Early', 'Mid', 'Late'))
metadata$PlantType <-
  factor(metadata$PlantType, levels = c(1, 2, 3, 4, 5) )

# this provides pathway abundance comparisons
fit_out <- maaslin3(input_data = pa.input,
                    input_metadata = metadata,
                    output = 'maaslin.pathway.planttype_and_stage.output',
                    formula = '~ PlantType + Stage + reads',
                    reference = 'Stage,Early',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1,
                    evaluate_only = "abundance",
                    warn_prevalence = FALSE,
                    coef_plot_vars = "Stage Mid,
                                      Stage Late",
                    heatmap_vars = "PlantType 2,
                                    PlantType 3,
                                    PlantType 4,
                                    PlantType 5"
                    )

# for both pathway and taxa,
# I used different combinations of formulas and variable plotting specifications to make multiple output folders
# the log file in each folder contains the exact formula used as well as values for the other variables
# taxa is the same for u50 or u90 runs. Pathway is different, so I ran those for both classification databases.

# this provides taxon abundance comparisons
fit_out.tax <- maaslin3(input_data = tax.input,
                        input_metadata = metadata,
                        output = 'maaslin.taxa.stage.output',
                        formula = '~ Stage + reads',
                        reference = 'Stage,Early',
                        normalization = 'TSS',
                        transform = 'LOG',
                        augment = TRUE,
                        standardize = TRUE,
                        max_significance = 0.1,
                        median_comparison_abundance = TRUE,
                        median_comparison_prevalence = FALSE,
                        max_pngs = 100,
                        cores = 1,
                        evaluate_only = "abundance",
                        warn_prevalence = FALSE,
                        coef_plot_vars = "Stage Mid,
                                          Stage Late"
                        )
```

For C4, I'm just copying in the project-specific code lines; most of it is standardized. The original script includes a way to change the column names of the metaphlan file to the short sample IDs, but I had already edited that file for the MaAsLin3 script. So that will be a convenient change for the future! I would prefer if I could show the names, especially for families, but the font has to be too small to fit them all, unfortunately.

```
setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann90/taxonomic-classifications")

clusters <- 7 # species
clusters <- 9 # genus; 7 probably would have been just as good
clusters <- 7 # family; anywhere from 6-9 would be fine
```

Again for C5 I'm copying in project-specific code lines. They are really handy to have if I need to recreate the code for a publication or a rerun.

```
metadata <- read.delim("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-metadata.txt", sep="\t")
metadata$PlantType <- as.character(metadata$PlantType)

setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-humann90/taxonomic-classifications")
```


## Module D - Co-Assembly of Sample Reads

Yet again I'm using the branch files here so I can correct or improve any of the code. This module will build a co-assembled metagenome for all samples in the project and back-align all samples to it.

### Script D1

The script creates the co-assembly. If the script times out before the assembly is complete, the code below can be rerun with the `--resume` flag. (It did time out, which makes sense with 15 samples!)

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D1.megahit-assemble.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly \
        --resume
```

### Script D2

The default assembly name given in step D1 is "final.contigs.fa"; it's also the default value in script D2, but I specified it just to be sure.

This script includes both the bowtie2 indexing and the creation of an assembly graph and reports with gfastats. Once the indexing is complete, D3 (back-alignment) and F1 (contig annotation) can be started without waiting for the assembly graph to be completed.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D2.bowtie-index.sh \
        --assemblyDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly \
        --prefix "mgx.coassembly" \
        --assemblyName "final.contigs.fa"
```

### Script D3

This is to back align all samples to the assembly, which makes abundance profiling and functional comparison possible downstream. G1 for binning can be started once the alignment is complete.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D3.bowtie-align.wrapper.sh \
        --fastqDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
        --assembly /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly \
        --bowtiePrefix /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/mgx.coassembly \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-alignment

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D3.bowtie-align.wrapper.sh \
        --fastqDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-fastq/cut-and-filtered/fastq \
        --assembly /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly \
        --bowtiePrefix /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/mgx.coassembly \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-alignment \
        --list "E3-cut"
```

## Module F - Annotation of Coassembly Contigs

### Script F1

This first script predicts protein sequences with CAT (which also provides taxonomic information), which I'll annotate with MicrobeAnnotator in the next step.

I noticed in the transplant side of things that this process hadn't completed all the way; there I could restart from alignment and not have to rerun F2, but here very very unfortunately even the predicted protein FAA and GFF files were incomplete. So I'm deleting it all and starting over. Then I'll have to start F2 over again also.

The predicted proteins and alignment completed but it needs to be run one more time to finish classification. I'm running into OOM kill errors with 320G for this classification! So I've edited the script to call a high memory node requesting 1024G; hopefully that helps.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F1.cat-contig-annotate.sh \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F1.cat-contig-annotate.sh \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation \
        --catProts /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/out.CAT.predicted_proteins.faa \
        --catAlign /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/out.CAT.alignment.diamond
```

### Script F2

This script assigns functional annotations to the proteins predicted by CAT using MicrobeAnnotator. I made edits in the code so I could run this on Phoenix, since Sol maintenance is coming up! Unfortunately it didn't work very well on Phoenix, so I'm restarting it on Sol. The files were actually still quite large after being split into 10 parts, so I split the original into five new files and I'm going to run the script on each part individually (so, 50 total files per the one fasta). Hopefully that makes this process complete sooner.

```
list=$(find ./ -name "*part*.0*.faa" | while read F; do basename $F; done)

for i in $list
do
echo $i
  sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.sh \
        /data/gencore/databases/microbe-annotator \
        /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/"$i" \
        "FALSE"
done
```

Unfortunately, 11 of the parts didn't complete within the four day time limit, so I'm needing to rerun them. I think there was some sort of error with the script hanging on the node they were on, because the process log file was never generated, so I can't even use the 'continue' option to speed things up.

```
list=$(find ./ -maxdepth 1 -name "*part*.0*.faa" | while read F; do basename $F; done) #(I moved the 11 files to a new dir to do this)

for i in $list
do
echo $i
  sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.sh \
        /data/gencore/databases/microbe-annotator \
        /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/"$i" \
        "FALSE"
done
```

Three of the files seemed to be hanging on the HMM profiling again, so I simultaneously am running them again into a different output folder.

```
list="out.CAT.predicted_proteins.part2.07.faa out.CAT.predicted_proteins.part3.04.faa out.CAT.predicted_proteins.part5.07.faa"

for i in $list
do
echo $i
  sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.sh \
        /data/gencore/databases/microbe-annotator \
        /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/"$i" \
        "FALSE"
done
```

It's finally complete! Yay!

### Scripts F3 and F4

I ran F3 manually to merge all 50 of the individual MicrobeAnnotator outputs. So now I can proceed with the remaining scripts in the analysis. For F4 I used the following input to the batch script:

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F4.ko_mapper.sh \
        -t /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        -s /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline \
        -c /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/final.contigs.fa \
        -o /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation
```

## Module G - Binning for MAGs

### Script G1

This script bins contigs from the coassembly together to predict metagenome-assembled genomes. It requires the completed alignments from D3 which is taking a while sigh.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G1.metabat-bin.sh \
        --bamDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-alignment \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bins
```

### Script G2

This script attempts to classify each bin using BAT, creating taxonomy predictions with scores.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G2.bat-bin-annotate.sh \
        --binDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bins/final.contigs.fa.metabat-bins-20240920_145447 \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bins
```

### Script G3

This script uses CheckM to classify each bin and summarize the quality of each MAG (homogeneity, contamination, and completeness). The classifications here are typically much more stringent than those from BAT, and the quantitative scores are good for selecting higher quality MAGs for polishing.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G3.checkm.sh \
        --binDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bins/final.contigs.fa.metabat-bins-20240920_145447 \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-bins
```

## Module H - Abundance Profiling

### Script H1

This script using the sample alignment files to the coassembly file and the CAT predicted protein GFF file to characterize the abundance of those proteins within the individual samples. I did have to run this twice due to an alignment error the first time.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_H/H1.abundance-profiling.wrapper.sh \
        --bamDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-alignment \
        --gffFile /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/out.CAT.predicted_proteins.gff \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-abundance-profiling
```

### Script H2

This script hasn't been streamlined for parameter inputs; the folder I used was `/data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-abundance-profiling`, the sample HTSeq counts file for header, etc. was `E1-cut_full.htseqcounts.txt`, and the script call was `sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_H/H2.merge-htseq.sh`. The goal of this script is to create a gene count matrix file for all genes present in each sample using the coassembly as a reference. Here is the exact code used for this project (after loading environments, SBATCH parameters, etc):

```
cd /data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-abundance-profiling

# change the sample prefix to any sample in the dataset
awk -F\\ '{print $1}' E1-cut_full.htseqcounts.txt > transcripts.txt
sed -i "1s/^/transcript_id \n/" transcripts.txt

for i in $(find ./ -type f -name "*htseqcounts.txt" | while read F; do basename $F; done | cut -d "." -f 1 )
do
  echo "$i"
  awk -F 't' '{print $NF}' "$i.htseqcounts.txt" > "$i.counts-only.txt"
  sed -i "1s/^/$i \n/" "$i.counts-only.txt"
done

paste transcripts.txt $(find ./ -type f -name "*.counts-only.txt" | sort | uniq ) > gene.count.matrix.tsv
```

### Script H4 Annotate HTSeq (lol even though it is the third H script)

I haven't run this yet; I'm waiting for MicrobeAnnotator to complete (running on Phoenix 10/11). gff is the original predicted proteins gff file from CAT (output of F1), matrix is the gene count matrix from HTseq abundance profiling (output of H2), anns is the annotated protein fasta file from microbe annotator (output of F3), and all the other variable files are generated within this script.

```
gff="/data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/out.CAT.predicted_proteins.gff"
newgff="nocomments.out.CAT.predicted_proteins.gff"
sortedgff="sorted.nocomments.out.CAT.predicted_proteins.gff"
annotations="suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"

anns="/data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-contig-annotation/all.out.CAT.predicted_proteins.faa.annot"
sortedanns="sorted.all.out.CAT.predicted_proteins.faa.annot"

matrix="/data/gencore/analysis_projects/7232686_Ball_WGS/carlini-transect/transect-abundance-profiling/gene.count.matrix.tsv"
tabbedGFF="gene.count.matrix_tabbedgff.tsv"
sortedTabbedGFF="sorted.gene.count.matrix_tabbedgff.tsv"

sortedAnnotations="sorted.suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"
geneInfo="complete-annotations-nocounts.tsv"

sed '/^#/ d' < "$gff" > "$newgff"

# sort the files
LC_ALL=C
LANG=C
sort -k 9 "$newgff" > "$sortedgff"
sort -k 1 "$anns" > "$sortedanns"

awk '{ split($9, a, ";"); print a[1], $1, $2, $3, $4, $5, $6, $7, $8, $9; }' "$sortedgff" \
| awk '{ split($1, a, "="); print a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10; }' \
| awk '{ split($1, a, "_"); print $1, $2 "_" a[2], $2, $3, $4, $5, $6, $7, $8, $9, $10; }' > "$annotations"

LC_ALL=C
LANG=C
join -o auto -a1 -e "null" <(sort "$matrix") <(sort "$annotations") | tr ' ' '\t' > "$tabbedGFF"

LC_ALL=C
LANG=C

sort -k 17 "$tabbedGFF" > "$sortedTabbedGFF"
join -o auto -a1 -e "null" -1 17 -2 1 -t $'\t' "$sortedTabbedGFF" "$sortedanns" > annotated.gene.count.matrix.tsv
# the value for sorting the tabbed gff file should be two more than the number of samples in the matrix

sort -k 2 "$annotations" | tr ' ' '\t' > "$sortedAnnotations"
join -o auto -a1 -e "null" -1 2 -2 1 -t $'\t' "$sortedAnnotations" "$sortedanns" > $geneInfo

# the geneInfo file is good for annotating paired MTX data!

# add header line to geneInfo file
sed -i '1i\ORFid\tGFFid\tContigID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tgffAttributes\tproteinID\tproteinDesc\tkeggID\tkeggDesc\ttaxonomy\tGO:MF\tGO:CC\tGO:BP\tIPRid\tpfamID\tECnum\tdatabase' $geneInfo
```
