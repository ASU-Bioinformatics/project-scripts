# Transect Metagenomic Workflow

## Module A

For this module, I'm starting with A1.trim-and-cut from the branch version of these scripts, to take advantage of cutadapt's superior adapter trimming while preserving as many reads as possible.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_A/A1.trim-and-cut.sh \
       /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq

mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/original
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-only
mkdir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-unpaired

mv *cut_SCT* /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-only
mv *cut_SQP* /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered
mv *cut_SUN* /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-unpaired
mv *fastq.gz /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/original
```

After adapter sequence and low quality reads are removed, we can move on to assembly-free analysis with Kraken, assembly-free analysis with HUMAnN, and coassembly. Since these are soil samples, we don't need to align to a reference genome to remove host reads.

## Module B - Assembly-Free Analysis with Kraken

I am also using the branch files for this module, because the original files don't have help messages and I'd like to add those and potentially clean up other aspects of the code as I go.

### Script B1

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B1.kraken.sh \
 --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
 --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-kraken \
 --dbDir /data/gencore/databases/kraken/k2_pluspf
```

### Script B2

Here, I'm choosing not to use the levels option because I'm satisfied with the default option (all taxonomic levels).

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B2.bracken.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-kraken/kraken-reports \
  --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bracken \
  --dbDir /data/gencore/databases/kraken/k2_pluspf
```

### Script B3

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B3.bracken-diversity.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bracken/bracken-raws \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bracken/bracken-diversity
```

### Script B4

I am using the default installations of KronaTools and KrakenTools to create the plots with this script, so they're not specified here in the command line call.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_B/B4.krona-plots.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bracken/bracken-reports \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-krona/krona-tables \
        --htmlDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-krona/krona-htmls
```

### Script B5 (in R)

I'm recording the custom parts of the script here in case I need to rerun anything! This code is specifically set up for species-level alpha diversity visualization, but I also ran the other taxonomic levels by replacing all instances of 'species' with the other taxonomic levels names respectively.

```
setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bracken/bracken-diversity/species")
div <- read.delim("species_alpha-diversity.txt")
meta <- read.delim("../../transplant-metadata.txt")
sampleids <- c("Bare", "Des", "Nos", "Pol", "San")
metadata <- c("PlantType")
# no numeric metadata set; the int values in PlantType are discrete and not a numeric series or continuous variables
```


## Module C - Assembly-Free Analysis with HUMAnN

Again I'm using the branch files here so I can correct or improve any of the code.

### Script C1

This script is the heart of the HUMAnN and MetaPhLAN pipeline. Because the fastq files haven't yet been concatenated, I'm using the default value for `--already-concatenated` which ensures that the forward and reverse reads are concatenated into a single file prior to classification. This script will be run twice - once for the uniref50 clustered database, and once for the uniref90 clustered database.

Because the database setting is global, DO NOT start running the u50 script until the u90 script is complete! If running multiple projects at once, DO NOT run u50 for any of them until all u90 runs are complete!

Additionally, a 'resume' option can be specified if a sample or samples aren't completed in the initial run through the pipeline; I had to do this for sample Des-cut and it is still taking a long time. I think it is because the file size is 3-4 times larger than any of the other samples!

I was able to get the analyses to complete more quickly by increasing the number of threads called by humann3 and the number of nodes requested from the server.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
  --dbType "u90" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann90

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq/concatenated \
  --already-concatenated --dbType "u90" --resume --list "Des-cut" \
  --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann90

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann50

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq/concatenated \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann50 \
  --already-concatenated --resume --list "Bare-cut Des-cut"

# this rerun I accidentally started on Phx. It should still work, just a reminder to check that server for job status.
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq/concatenated \
  --dbType "u50" --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann50 \
  --already-concatenated --resume --list "Des-cut"
```

### Other C Scripts (2-5)
After this script completes, I run the code in `C2.humann-manipulations.sh` interactively to standardize file formatting and location within the output directory.

Scripts `C3.maaslin.R`, `C4.metaphlan-heatmap.R`, and `C5.metaphlan-alpha-diversity.R` are run interactively in a local R session, allowing figures to be fine-tuned and giving me a look at what's going on with the data.

Since a new version of MaAsLin is available now, I added a new script, `C3A.maaslin3.R` to Module_C and ran it. It allowed me to examine compositional abundance and prevalence comparisons for both pathways and taxa, which is really nice! A lot of custom variables are used in this script, from the input data to the model fit formula, so I'm including the entire R script here.

The same script was used for u50 as well, just changing the appropriate input/output folder names.

```
for (lib in c('maaslin3', 'dplyr', 'ggplot2', 'knitr', 'kableExtra')) {
  suppressPackageStartupMessages(require(lib, character.only = TRUE))
}

#### DEFINE CUSTOM VARIABLES ####

# this is only for nonpaired data. Use MTX model scripts for paired metatranscriptomic data (section J)
# to prepare abundance file for R, remove hash from Pathway column and shorten the sample names

metaDir <- "/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/"
metaFile <- paste0(metaDir,"transplant-metadata.txt")
projectDir <- "/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann90/"
outputDir <- paste0(projectDir, "pathways/differential-pathways")
inputDir <- paste0(projectDir, "pathways/pathway-abundance/")
taxaDir <- paste0(projectDir, "taxonomic-classifications/")
normalizedAbundance <- paste0(inputDir, "normalized_pathabundance_unstratified.tsv")
metaphlanTaxonomy <- paste0(taxaDir, "merged_metaphlan_table_species.txt")

#### SET DIRECTORY AND UPLOAD METADATA ####

setwd(outputDir)

metadata <- data.frame(read.delim(metaFile, sep="\t"))

rownames(metadata) <- metadata[,1]
metadata$reads <- c(40442320, 100208623, 28585644, 25908193, 27685642)
metadata <- metadata[,-1]

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
metadata$PlantType <-
  factor(metadata$PlantType, levels = c("Bare", "Des", "Nos", "Pol", "San") )

# this provides pathway abundance comparisons
fit_out <- maaslin3(input_data = pa.input,
                    input_metadata = metadata,
                    output = 'u90.maaslin.output',
                    formula = '~ PlantType + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)

# this provides taxon abundance comparisons
fit_out.tax <- maaslin3(input_data = tax.input,
                    input_metadata = metadata,
                    output = 'u90.maaslin.taxon.output',
                    formula = '~ PlantType + reads',
                    normalization = 'TSS',
                    transform = 'LOG',
                    augment = TRUE,
                    standardize = TRUE,
                    max_significance = 0.1,
                    median_comparison_abundance = TRUE,
                    median_comparison_prevalence = FALSE,
                    max_pngs = 100,
                    cores = 1)
```

For C4, I'm just copying in the project-specific code lines; most of it is standardized. The original script includes a way to change the column names of the metaphlan file to the short sample IDs, but I had already edited that file for the MaAsLin3 script. So that will be a convenient change for the future! I would prefer if I could show the names, especially for families, but the font has to be too small to fit them all, unfortunately.

```
rowMD <- c("blue", "red", "green", "orange", "purple")
setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann90/taxonomic-classifications")

clusters <- 3 # species
clusters <- 5 # genus; 3 probably would have been just as good
clusters <- 4 # family; anywhere from 3-6 would be fine
```

Again for C5 I'm copying in project-specific code lines. They are really handy to have if I need to recreate the code for a publication or a rerun. With only one sample per sample group, I can't graph alpha diversity for this project.

```
metadata <- read.delim("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-metadata.txt", sep="\t")

setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-humann90/taxonomic-classifications")
```


## Module D - Co-Assembly of Sample Reads

Yet again I'm using the branch files here so I can correct or improve any of the code. This module will build a co-assembled metagenome for all samples in the project and back-align all samples to it.

### Script D1

The script creates the co-assembly. If the script times out before the assembly is complete, the code below can be rerun with the `--resume` flag. With only five samples, this assembly was completed within the original time allotment.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D1.megahit-assemble.sh \
        --inputDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly
```

### Script D2

The default assembly name given in step D1 is "final.contigs.fa"; it's also the default value in script D2, but I specified it just to be sure.

This script includes both the bowtie2 indexing and the creation of an assembly graph and reports with gfastats. Once the indexing is complete, D3 (back-alignment) and F1 (contig annotation) can be started without waiting for the assembly graph to be completed.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D2.bowtie-index.sh \
        --assemblyDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly \
        --prefix "mgx.coassembly" \
        --assemblyName "final.contigs.fa"
```

### Script D3

This is to back align all samples to the assembly, which makes abundance profiling and functional comparison possible downstream. It didn't to work for Des-cut (always the problem sample) because I inadvertently moved the fastq files instead of copying them to try a fresh starts for the humann90 run, so I'm rerunning it now. G1 for binning can be started once the alignment is complete.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D3.bowtie-align.wrapper.sh \
        --fastqDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
        --assembly /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly \
        --bowtiePrefix /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/mgx.coassembly \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-alignment

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_D/D3.bowtie-align.wrapper.sh \
        --fastqDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-fastq/cut-and-filtered/fastq \
        --assembly /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly \
        --bowtiePrefix /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/mgx.coassembly \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-alignment \
        --list "Des-cut"
```

## Module F - Annotation of Coassembly Contigs

### Script F1

This script predicts ORFs, aligns them to the CAT database, and then predicts taxonomic data from the alignments.

The first time I ran this I thought it was complete because the predicted proteins and GFF file had been created. But the alignment and classification files hadn't been generated due to a timeout, unfortunately. I only realized when script F4 had errors. So I'm rerunning now with a longer time limit from the alignment process forward.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F1.cat-contig-annotate.sh \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F1.cat-contig-annotate.sh \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation \
        --catProts /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.faa
```

### Script F2

This script assigns functional annotations to the proteins predicted by CAT using MicrobeAnnotator. Chunk 06 didn't complete annotation so I've resumed that analysis.

Then, there was an error on the resume for chunk six so maybe I need to set this up differently hmm... I'm trying to just get into the annotator script directly instead of using the wrapped; had to manually update the outDir to do so, so I'll need to reset that for the next batch after.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.wrapper.sh \
        --databaseDir /data/gencore/databases/microbe-annotator \
        --proteinFasta /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.faa --split

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.wrapper.sh \
        --databaseDir /data/gencore/databases/microbe-annotator \
        --proteinFasta /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.faa --continue --list "6"

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F2.microbe-annotator.sh \
       /data/gencore/databases/microbe-annotator \
       /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.06.faa \
       "FALSE"
```

### Script F3 and F4

F2 took forever but finally fragment 6 completed annotations!

*F4 needs to be rerun once F1 completes its rerun* - and apparently the rerun timed out and I just caught that. So here we go again.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F3.microbe-annotator-merge.sh \
        /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F4.ko_mapper.sh \
        -t /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        -s /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline \
        -c /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/final.contigs.fa \
        -o /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_F/F4.ko_mapper.sh \
        -t /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        -s /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline \
        -c /data/gencore/analysis_projects/7232686_Ball_WGS/transplant-to-return/4.metagenome-assembly/final.contigs.fa \
        -o /data/gencore/analysis_projects/7232686_Ball_WGS/transplant-to-return/5.assembly-annotation
```

## Module G - Binning (and Eventually MAG Reporting)

### Script G1

This script outputs the bins in the alignment folder; I typically just move them to a new output directory manually. Definitely too tired now to improve the code lol.

Actually I have fixed this now so an output folder can be specified! Much nicer. It does create a subfolder for the bins themselves (the subfolder and a depth txt file are put into the specified output directory from the script).

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G1.metabat-bin.sh \
        --bamDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-alignment \
        --contigs /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-assembly/final.contigs.fa \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bins
```

### Script G2

This script uses BAT (the sister of CAT) to annotate and classify the bins. The database and taxonomy directories are the same as those used for F1 CAT annotation of the coassembly contigs.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G2.bat-bin-annotate.sh \
        --binDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bins/final.contigs.fa.metabat-bins-20240905_204215 \
        --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
        --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bins
```

### Script G3

This script uses CheckM to classify each bin and summarize the quality of each MAG (homogeneity, contamination, and completeness). The classifications here are typically much more stringent than those from BAT, and the quantitative scores are good for selecting higher quality MAGs for polishing.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_G/G3.checkm.sh \
        --binDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bins/final.contigs.fa.metabat-bins-20240905_204215 \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-bins
```

## Module H - Contig-Based Abundance Profiling

### Script H1

Using the contig annotations from CAT and the sample alignments to the contig assembly, this script counts the number of reads corresponding to each predicted protein for comparative analysis downstream.

```
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway-release1.0.0-branch/Module_H/H1.abundance-profiling.wrapper.sh \
        --bamDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-alignment \
        --gffFile /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.gff \
        --outDir /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-counts
```

### Script H2

This script isn't set up for input parameters - the code itself has to be altered. So I'm copying out that sbatch code here for future reference.

```
module load mamba/latest

source deactivate
source activate /data/biocore/programs/mamba-envs/htseq-env

# move all htseqcounts.txt files to a new directory before running
# then cd to that directory

cd /data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-counts

# change the sample prefix to any sample in the dataset
awk -F\\ '{print $1}' San-cut_full.htseqcounts.txt > transcripts.txt
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

```
# gff is the original predicted proteins gff file from CAT (output of F1)
# matrix is the gene count matrix from HTseq abundance profiling (output of H2)
# anns is the annotated protein fasta file from microbe annotator (output of F3)
# all the other variable files are generated within this script.

gff="/data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/out.CAT.predicted_proteins.gff"
newgff="nocomments.out.CAT.predicted_proteins.gff"
sortedgff="sorted.nocomments.out.CAT.predicted_proteins.gff"
annotations="suffix.sorted.nocomments.out.CAT.predicted_proteins.gff"

anns="/data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-contig-annotation/all.out.CAT.predicted_proteins.faa.annot"
sortedanns="sorted.all.out.CAT.predicted_proteins.faa.annot"

matrix="/data/gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-counts/gene.count.matrix.tsv"
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
#cat $matrix | tr '\t' ' ' > spaced.matrix.tsv
join -a1 -o auto -e "null" <(sort "$matrix") <(sort "$annotations") | tr ' ' '\t' > "$tabbedGFF"

LC_ALL=C
LANG=C

sort -k 7 "$tabbedGFF" > "$sortedTabbedGFF"
join -o auto -a1 -e "null" -1 7 -2 1 -t $'\t' "$sortedTabbedGFF" "$sortedanns" > annotated.gene.count.matrix.tsv
# the value for sorting the tabbed gff file should be two more than the number of samples in the matrix

sort -k 2 "$annotations" | tr ' ' '\t' > "$sortedAnnotations"
join -o auto -a1 -e "null" -1 2 -2 1 -t $'\t' "$sortedAnnotations" "$sortedanns" > $geneInfo

# add header line to geneInfo file
sed -i '1i\ORFid\tGFFid\tContigID\tsource\ttype\tstart\tend\tscore\tstrand\tphase\tgffAttributes\tproteinID\tproteinDesc\tkeggID\tkeggDesc\ttaxonomy\tGO:MF\tGO:CC\tGO:BP\tIPRid\tpfamID\tECnum\tdatabase' $geneInfo

# after annotating, use R to create a normalized matrix
# R code: df[-1] <- sapply(df[-1], prop.table) * 100
```

### Script H4 Extract Keggs

Again this is a manual script, so I'm including the exact code we used here.

```
sed -i.txt '1i contigID\ttranscriptID\tBare\tDes\tNos\tPol\tSan\tseqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute\tncbiID\tncbiName\tkeggID\tkeggName\ttaxonomy\tgoMF\tgoCC\tgoBP\tinterproscanID\tpfamID\tuniProtEnzyme\tdatabase' annotated.gene.count.matrix.tsv

declare -A sampleArray

sampleArray+=( ["3"]=Bare
               ["4"]=Des
               ["5"]=Nos
               ["6"]=Pol
               ["7"]=San )

# the column number is 14 + the sample number
for key in ${!sampleArray[@]}; do
    echo ${key} ${sampleArray[${key}]}
    awk -v col="${key}" -F $'\t' '($col!=0 && $19!="null" && $19!="NA") {print $19}' annotated.gene.count.matrix.tsv > "${sampleArray[${key}]}".keggs.tsv
done

kolist=$(find ./ -maxdepth 1 -type f -name "*.keggs.tsv" | while read F; do basename $F; done)

python /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline/ko_mapper.py \
  -i $kolist \
  -p ko_map --cluster rows

# find the number of unique Kegg terms identified in at least one sample (the column number is 14 + the number of samples)
awk -F '\t' '{print $19}' annotated.gene.count.matrix.tsv | grep -o -E 'K[[:digit:]]+' | sort | uniq | wc -l
#ball transplant = 9300

# find the number of unique GO MF terms identified in at least one sample (the column number is 17 + the number of samples)
awk -F '\t' '{print $22}' annotated.gene.count.matrix.tsv | grep 'GO' | tr ' ' '\n' | sort | uniq | wc -l
# ball transplant = 4113

# find the number of unique GO CC terms identified in at least one sample (the column number is 18 + the number of samples)
awk -F '\t' '{print $23}' annotated.gene.count.matrix.tsv | grep 'GO' | tr ' ' '\n' | sort | uniq | wc -l
# ball transplant = 1336

# find the number of unique GO BP terms identified in at least one sample (the column number is 19 + the number of samples)
awk -F '\t' '{print $24}' annotated.gene.count.matrix.tsv | grep 'GO' | tr ' ' '\n' | sort | uniq | wc -l
# ball transplant = 6540

# find the number of unique taxa identified in at least one sample (the colum number is 16 + the number of samples)
awk -F '\t' '{print $21}' annotated.gene.count.matrix.tsv | sort | uniq | wc -l
# ball transplant = 33,323
```

The genetic potential for this dataset includes 9300 unique Kegg terms, 4113 unique GO molecular function terms, 1336 unique GO cell component terms, and 6540 unique GO biological process terms. 33,323 unique taxa were identified via representative proteins. If I don't already have this in R with the normalized annotated matrix, I want to get the relative proportions of those taxa for each sample as I think it would be very valuable. (I do not have this; I would need to write a script for it. Maybe tomorrow? The annotated gene count matrix may also be too large for my local R setup to handle. Script H5 only works with the Kegg module table, which is a lot smaller. So I would probably need to do it through a Sol R interactive session.)

### Script H5: Kegg Heatmap in R

As with most of my R scripts, this is a lot more hands-on, so I'm copying the script as I ran it here in case it's needed for future reference.

```
library("gplots")
library(tidyverse)
library(vegan)
library(ggstatsplot)

setwd("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-counts")

keggs <- read.delim("ko_map_module_completeness.tab", sep="\t")

kegg.matrix <- keggs[,c(2:8)]
rownames(kegg.matrix) <- kegg.matrix[,1]
kegg.matrix <- kegg.matrix[,2:7]
kegg.matrix$variance <- apply(kegg.matrix[,2:6],1,var)
kegg.matrix <- kegg.matrix[order(kegg.matrix$variance, decreasing = TRUE),]
kegg.modules <- kegg.matrix[1:100,]$pathway.group
kegg.matrix <- as.matrix(kegg.matrix[1:100,2:6])


clustRowBar <- viridis::turbo(length(unique(kegg.modules)), begin=0, end=1)

kegg.modules <- as.data.frame(kegg.modules)
kegg.modules$kegg.modules <- as.factor(kegg.modules$kegg.modules)
kegg.modules$kegg.moduleNames <- as.factor(kegg.modules$kegg.modules)
levels(kegg.modules$kegg.modules) <- as.factor(clustRowBar)
clustRowBar <- as.vector(kegg.modules$kegg.modules)
labels <- c("San", "Des", "Bare", "Nos", "Pol")

pdf("kegg.heatmap.pdf", height=24, width=20)
par(mar = c(2, 2, 16, 2),                                  # Specify par parameters
    xpd = TRUE)
heatmap.2(kegg.matrix[1:100,1:5],
          Rowv=TRUE,
          Colv=TRUE,
          col=redgreen(100),
          scale="col",
          dendrogram = "col",
          #margins = c(1, 1),
          keysize=0.6,
          cexCol = 1,
          cexRow = 0.9,
          labRow = rownames(kegg.matrix[1:100,1:5]),
          labCol = labels,
          main = "100 Most Variable Kegg Modules by Completeness",
          RowSideColors = clustRowBar,
          trace = "none",
          margins=c(6,44))
#legend(x="top", inset = c(-0.06, -0.12), cex=0.9,
 #      legend = levels(kegg.modules$kegg.moduleNames),
  #     fill = clustRowBar, ncol = 3)
dev.off()

#### TMM Normalization of Gene Matrix ####
library(edgeR)
y<-DGEList(counts=kegg.matrix)
y<-calcNormFactors(y)
d<-cpm(y) # since calcNormFactors was run first, this is TMM rather than just counts-per-million output
write.table(d,"TMM.normalized.counts.txt",quote=FALSE,sep="\t",row.names=TRUE)

#### PCA PLOT BY GENE COUNT MATRIX ####
metadata <- read.delim("/Volumes/Gencore/analysis_projects/7232686_Ball_WGS/transplant/transplant-metadata.txt", sep="\t")
metadata <- metadata[order(metadata$SampleID),]

genecount_df <- d
genecount_df <- genecount_df[,order(colnames(genecount_df))]

#colnames(genecount_df) <- c(metadata$Sample.ID)

rownames(metadata) <- metadata$Sample.ID
metadata <- as.data.frame(metadata[,-1])
colnames(metadata) <- "PlantType"

gene_mat = genecount_df |>
  #column_to_rownames("transcript_id") |>
  as.matrix() |>
  t()

dist_mat = vegdist(gene_mat)

cmd_res = cmdscale(dist_mat,
                   k = (nrow(gene_mat) - 1),
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


wa_data = wascores(cmd_res$points[,1:2], gene_mat) |>
  as_tibble(rownames = 'geneID')

wa_data
```
