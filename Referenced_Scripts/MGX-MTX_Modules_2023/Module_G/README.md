# Module G
This module uses [METABAT](https://bitbucket.org/berkeleylab/metabat/src/master/,"Title") to bin contigs into potential metagenome-assembled genomes (MAGs), using the alignment data from each sample against the contig coassembly as a guide. CheckM and BAT are used to classify each bin.

Two optional scripts allow for additional evaluation of high-scoring bins if a goal is to report MAGs or submit bins as MAGs on NCBI. These scripts do require the additional work of identifying and downloading the appropriate comparative genomes (other strains in the same species, or species in the same genus) to use.

This module is intended primarily for MGX data; due to its utilization of average base coverage in the alignments to separate contigs into bins, METABAT2 is not expected to perform well for MTX data.

## G1.metabat-bin.sh
METABAT2 is an unsupervised algorithm to cluster contigs based on probabilistic distance calculated from tetra-nucleotide frequency and average base coverage in a sample. Because it is not limited by reference genomes, it has the potential to capture environmental genomes that are previously uncharacterized.

To determine and use the average base coverage, METABAT2 requires sorted BAM alignment files (of DNA reads to the contig assembly file) as input.

A sample call:

    sbatch ./G1.metabat-bin.sh \
      --bamDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/metabat-bins \
      --contigs /Volumes/Gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa

Input parameters:
* `-b,--bamDir`: the directory containing the sorted bam files for sample alignments to the contig coassembly. The folder containing the assembled bins will also be output here.
* `-c,--contigs`: the contig coassembly fasta file.

## G2.bat-bin-annotate.sh
As with Module F for contig annotations, we use the [CAT/BAT](https://github.com/dutilh/CAT,"Title") tool for protein predictions and taxonomic classifications within bins. The tool relies on the pre-prepared CAT database and taxonomy as a reference.

On Sol, this script requires the Conda environment `/data/biocore/programs/mamba-envs/catbat-env` installed on the `/data/biocore` server. If this is not an option, an equivalent Conda environment can be installed using the included .yml file.

A sample call:

    sbatch ./G2.bat-bin-annotate.sh \
      --binDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins/metabat-bins \
      --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
      --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins

Input parameters:
* `-b,--binDir`: the path to the directory containing bins in fasta format
* `-d,--databaseDir`: the path to the directory containing the CAT database files
* `-t,--taxDir`: the path to the directory containing the CAT taxonomy files
* `-o,--outDir`: the path to the directory where the output files should be written

This script produces a log file for the run, a text file containing the top hit for each ORF within each bin; a [Diamond](https://github.com/bbuchfink/diamond,"Title") alignment file for the bins against the database; nucleotide fasta, amino acid fasta, and gff files for the concatened set of predicted proteins; and classification tables with numeric taxa IDs or human readable taxa names.

## G3.checkm.sh
[CheckM](https://github.com/Ecogenomics/CheckM/wiki,"Title") evaluates the quality of each bin as a unique genome, estimating the completeness and contamination for each bin. Additionally, it predicts taxonomy for the bin. CheckM has much more conservative taxonomic assignments than BAT, tending to assign bins to higher taxonomic levels. However, the completeness and contamination estimates make it an invaluable tool.

On Sol, this script can be run from within the Conda environment `/data/biocore/programs/mamba-envs/checkm-env` on the core's data server. However, if a local installation is needed the included yml file can be used to emulate the original environment.

A sample call:

    sbatch ./G3.checkm.sh \
      --binDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins/metabat-bins \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/coassembly-bins

Input parameters:
* `-b,--binDir`: the path to the directory containing bins in fasta format
* `-o,--outDir`: the path to the directory where the output files should be written

Output files include an Excel file containing the classifications and quality scores for each bin, as well as intermediate files produced during the run and a log file for the process.
