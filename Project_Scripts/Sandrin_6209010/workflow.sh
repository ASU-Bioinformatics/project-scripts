# commands run for Keaton's MGX/MTX samples, beginning with raw short-read fastq files.

# Module A, script A1; no sbatch command, pathnames were updated manually

# Module A, script A2, remove human contamination from DNA samples
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A2.remove-host.wrapper.sh \
  --assembly /data/gencore/databases/reference_genomes/human/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome \
  --fqDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/fastq/adapter-trimmed/fastq/ \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
  --use "DIRECTORY"

# Module A, script A3, remove human contamination from RNA samples
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_A/A3.remove-host-rna.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/fastq/adapter-trimmed/fastq \
  --refDir /data/gencore/databases/reference_genomes/human/Ensembl_GRCh38_102_agaveSTAR2.7.3/star_2.7.10a_indexes \
  --tempDir /scratch/kawoodbu \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq \
  --use "DIRECTORY"

# Module B, script B1, classify DNA reads with Kraken following removal of host reads
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_B/B1.kraken.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/kraken-outputs \
  --dbDir /data/gencore/databases/kraken/k2_pluspf

# Module B, script B2, estimate abundance at each taxonomic level with Bracken
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_B/B2.bracken.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/kraken-outputs/kraken-reports \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs \
  --dbDir /data/gencore/databases/kraken/k2_pluspf

# Module B, script B3, calculate alpha and beta diversity from Bracken abundance estimates
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_B/B3.bracken-diversity.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-raws \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-diversity

# Module B, script B4, create Krona HTML sunburst plots and accompanying tables
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_B/B4.krona-plots.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/bracken-outputs/bracken-reports \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/krona-outputs \
  --htmlDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/kraken-pipeline/krona-htmls \
  --krakenBin /data/biocore/programs/KrakenTools-1.2 \
  --kronaBin /data/biocore/programs/KronaTools-2.8.1/bin

# Module C, script C1, run basic humann and metaphlan analysis on DNA samples for UniRef90
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann90_pipeline \
  --dbType u90

# Module C, script C1, run basic humann analysis on RNA samples with DNA population metrics, UniRef90
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq/ \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann90-paired_pipeline \
  --dnaRef /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann90_pipeline \
  --paired --dbType u90

# Module C, script C2, normalize, merge, and unstratify humann/metaphlan output for DNA samples, UniRef90

# Module C, script C2, normalize, merge, and unstratify humann/metaphlan output for RNA samples, UniRef90

# Module C, script C4, create heatmap for metaphlan taxonomic data (not run for paired RNA, same for UniRef50 and 90)

# Module C, script C5, plot PCA plots and alpha diversity box-violin plots for metaphlan data (only needed once in the project)

# Module C, script C1, run basic humann and metaphlan analysis on DNA samples for UniRef50
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq/concatenated \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann50_pipeline \
  --dbType u50 --concat

# Module C, script C1, run basic humann analysis on RNA samples with DNA population metrics, UniRef50
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_C/C1.humann-metaphlan.wrapper.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq/concatenated \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/humann50-paired_pipeline \
  --dnaRef /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/humann50_pipeline \
  --concat --paired --dbType u50

# Module C, script C2, normalize, merge, and unstratify humann/metaphlan output for DNA samples, UniRef50

# Module C, script C2, normalize, merge, and unstratify humann/metaphlan output for DNA samples, UniRef50

# Module D, script D1, assemble reads into contigs with MEGAHIT
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D1.megahit-assemble.sh \
  --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly

# Module D, script D2, build Bowtie2 Index for contig coassembly
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D2.bowtie-index.sh \
  --assemblyDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly \
  --assemblyName final.contigs.fa --prefix coassembly

# Module I, script I1, build HiSat2 Index for contig coassembly, for RNA alignment
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I1.hisat-index.sh \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-index \
  --prefix final.contigs.reference \
  --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa

# Module I, script I2, align transcripts to coassembly using HiSat2 index
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I2.align-transcripts.wrapper.sh \
  --transcriptDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/fastq \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-alignments \
  --hisatIndex /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-index/final.contigs.reference

# Module F, script F1, predict proteins on contigs in coassembly with CAT
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F1.cat-contig-annotate.sh \
  --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
  --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
  --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

# Module I, script I3, count reads for each sample on the coassembly gff and alignments
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I3.abundance-profiling.wrapper.sh \
  --alignmentDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/hisat2-alignments \
  --gffFile /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/out.CAT.predicted_proteins.gff \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts

# Module F, script F2, annotate proteins predicted by CAT using MicrobeAnnotator
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F2.microbe-annotator.wrapper.sh \
  --databaseDir /data/gencore/databases/microbe-annotator \
  --proteinFasta /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/out.CAT.predicted_proteins.faa

# Module F, script F3, merge annotations
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F3.microbe-annotator-merge.sh /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

# Module F, script F4, auxiliary scripts for creating some helpful files
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F4.ko_mapper.sh \
  --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
  --scriptDir /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline \
  --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

# Module D, script D3, sample alignment of DNA reads to bowtie2 index
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_D/D3.bowtie-align.wrapper.sh \
  --fastqDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
  --assembly /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/coassembly \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/bowtie2-alignments

# Module G, script G1, metabat binning given coassembly and back alignments
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_G/G1.metabat-bin.sh \
  --bamDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/bowtie2-alignments \
  --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa

# Module H, script H1, abundance profiling of sorted bam files with HTseq
bash /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_H/H1.abundance-profiling.wrapper.sh \
  --bamDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/bowtie2-alignments \
  --gffFile /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/out.CAT.predicted_proteins.gff \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/htseq-counts

# Module G, script G2, annotation and classification with BAT
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_G/G2.bat-bin-annotate.sh \
  --binDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins/metabat-bins \
  --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
  --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins

# Module G, script G3, classification with CheckM
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_G/G3.checkm.sh \
  --binDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins/metabat-bins \
  --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins

# Module I, script I4, merge htseq counts
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_I/I4.merge-htseq.sh \
  /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts

# Module K, script K1, subscript DESeq2 analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script deseq2 \
  --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression-redo \
  --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/gene.count.matrix.tsv \
  --comparisons /data/gencore/sftp/t_sandrin/6209010_Metaomics/7.contig-annotation/3.RNA-abundance-profiling-will-need-updating/differential-expression-will-need-updating/comparisons.csv

# Module K, script K1, subscript edgeR analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script edger \
  --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression-redo \
  --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/gene.count.matrix.tsv \
  --comparisons /data/gencore/sftp/t_sandrin/6209010_Metaomics/7.contig-annotation/3.RNA-abundance-profiling-will-need-updating/differential-expression-will-need-updating/comparisons.csv

# Module K, script K1, subscript NOIseq analysis
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_K/K1.rscripts.wrapper.sh \
  --script NOIseq \
  --directory /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/differential-expression-redo \
  --geneMatrix /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed/htseq-counts/gene.count.matrix.tsv \
  --comparisons /data/gencore/sftp/t_sandrin/6209010_Metaomics/7.contig-annotation/3.RNA-abundance-profiling-will-need-updating/differential-expression-will-need-updating/comparisons.csv

# Module J, scripts J1, J2, J3: paired normalization and metatranscriptomic evaluation of taxa/functional pathway differences

# non-modular scripts, including 16S scripts:
# 16S-qiime2-pt1.sh - denoising, feature identification, and basic sample stats for 16S data using Qiime2
# 16S-qiime2-pt2.sh - diversity stats, taxonomic classification and analysis, and ANCOM-BC using Qiime2
# 16S-presence-absence-reduced.R - dissimilarity indices and visualization for 16S data
# maaslin3.R - replaces maaslin2 scripts in referenced MGX script modules
