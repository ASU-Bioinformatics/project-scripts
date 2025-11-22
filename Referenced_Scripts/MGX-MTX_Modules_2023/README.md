# Standalone Metagenomics Pathway

This pathway is equipped for metagenomic, metatranscriptomic, paired metagenomic/metatranscriptomic, and viral enriched metagenomics and metatranscriptomic data. In the rest of the manual, MGX and MTX will be used for metagenomic and metatranscriptomic respectively.

## MGX Data Modules

To analyze a data set containing only MGX reads, you'll need the following script modules:
1. Module A: Quality filter and trim adapters
2. Module B: Kraken k-mer based taxonomic classification and ecological diversity
3. Module C: HUMAnN/MetaPhLAN taxonomic classification and functional analysis
4. Module D: Coassembly of reads in to contigs and back alignment of each sample to the coassembly
5. Module F: Annotation of contigs with functional analysis
6. Module G: Binning of contigs into potential MAGs, with optional phylogenetic comparison of individual MAGs to reference genomes
7. Module H: Abundance profiling of contig annotations across individual samples with functional analysis

## MTX Data Modules

To analyze a data set containing only MTX reads, you'll need the following script modules:
1. Module A: Quality filter and trim adapters
2. Module C: HUMAnN/MetaPhLAN taxonomic classification and functional analysis
3. Module D: Coassembly of reads in to contigs
4. Module F: Annotation of contigs with functional analysis
6. Module I: Alignment of transcripts to coassembly with abundance profiling
7. Module K: Differential expression analysis between sample groups

## Paired MGX/MTX Data Modules

To analyze a data set containing paired MGX and MTX reads, you'll need the following script modules:
1. Module A: Quality filter and trim adapters
2. Module B: Kraken k-mer based taxonomic classification and ecological diversity
3. Module C: HUMAnN/MetaPhLAN taxonomic classification and functional analysis
4. Module D: Coassembly of reads in to contigs and back alignment of each sample to the coassembly
5. Module F: Annotation of contigs with functional analysis
6. Module G: Binning of contigs into potential MAGs, with optional phylogenetic comparison of individual MAGs to reference genomes
7. Module H: Abundance profiling of contig annotations across individual samples with functional analysis
8. Module I: Alignment of transcripts to DNA coassembly with abundance profiling
9. Module J: MaAsLin multivariate modeling of transcript expression levels from MTX by population characteristics by MGS, comparing joint functional and taxonomic data between sample groups
10. Module K: Differential expression analysis between sample groups

## Viral-Enriched MGX Data Modules

To analyze a data set containing MGX reads enriched for viral sequences, you'll need the following script modules:
1. Module A: Quality filter and trim adapters
2. Module B: Kraken k-mer based taxonomic classification and ecological diversity
3. Module C: HUMAnN/MetaPhLAN taxonomic classification and functional analysis
4. Module D: Coassembly of reads in to contigs and back alignment of each sample to the coassembly
5. Module F: Annotation of contigs with functional analysis
6. Module G: Binning of contigs into potential MAGs, with optional phylogenetic comparison of individual MAGs to reference genomes
7. Module H: Abundance profiling of contig annotations across individual samples with functional analysis
8. Module E: Identification and classification of viral contigs and raw reads, with comparative diversity analysis

## Viral-Enriched MTX Data Modules

To analyze a data set containing MTX reads enriched for viral sequences, you'll need the following script modules:
1. Module A: Quality filter and trim adapters
2. Module C: HUMAnN/MetaPhLAN taxonomic classification and functional analysis
3. Module D: Coassembly of reads into contigs and back alignment of each sample to the coassembly
4. Module F: Annotation of contigs with functional analysis
6. Module I: Alignment of transcripts to coassembly with abundance profiling
8. Module K: Differential expression analysis between sample groups
7. Module E: Identification and classification of viral contigs and raw reads, with comparative diversity analysis
