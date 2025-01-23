# Module E: Viral Identification for Enriched Samples

Currently, this module is usable for MGX and MTX samples; it hasn't been expanded for paired MTX data. It approaches viral identification in several ways. First, [VIBRANT](https://github.com/AnantharamanLab/VIBRANT/tree/master,"Title") is used to identify viral contigs in the contig assembly file. Then, [vironomy](https://github.com/b-tierney/vironomy,"Title") is used to classify as many of those contigs as possible. Unclassified contigs might contain viral sequences integrated into other genomes, or viral sequences from uncharacterized viruses. Third, as a read-based vs. a contig-based strategy, all raw reads are aligned against the NCBI RefSeq viral database with [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml,"Title").

## E1.vibrant.sh

Vibrant takes the coassembled contig file as input and uses neural network machine learning based on protein signatures to predict viral contigs. Additionally, it analyzes the viral community function via identification of metabolic genes and pathways (focusing on KEGG terms).

A sample code:

    sbatch ./E1.vibrant.sh \
      --contigs /data/gencore/analysis_projects/6078853_Otak_RNA/megahit_alignments/final.contigs.fa \
      --databaseDir /data/gencore/databases/vibrant2/databases \
      --modelDir /data/gencore/databases/vibrant2/files \
      --outDir /data/gencore/analysis_projects/6078853_Otak_RNA/vibrant

Input parameters:
* `-c,--contigs`: the pathname for coassembly contig fasta file
* `-d,--databaseDir`: the directory in which the VIBRANT database is stored
* `-m,--modelDir`: the directory in which the VIBRANT model files are stored
* `-o,--outDir`: the directory to which VIBRANT output files should be written

Output:
There are a *lot* of output files from a VIBRANT run. This is a brief summary of the folders and several of the most important files.
* `VIBRANT_phages_<input_file>`: this folder contains FASTA, GenBank, and tables for identified viral contigs
* `VIBRANT_results_<input_file>`: this folder contains annotations, metabolic identifications, and function summaries for the data
    * `VIBRANT_complete_circular_<input_file>`: this file contains a list of all contigs that were determine to contain complete circular viral genomes
    * `VIBRANT_annotations_<input_file>`: this file contains a list of all identified contigs with KEGG and pfam annotations where present
    * `VIBRANT_AMG...`: these files contain count data for each identified auxiliary metabolic gene in the data set, the identified auxiliary metabolic genes on each contig, and the count data for each identified auxiliary metabolic pathway in the data set. Given that this is the characterization of the entire coassembly file, the file `VIBRANT_AMG_individuals_final.contigs.tsv` is most useful and can be combined with abundance profiling to examine the AMG profiles of each sample in the data set individually.
* `VIBRANT_HMM_tables_<type>_<input_file>`: the parsed and unformatted folders contain the raw HMM tables used for classification and functional analysis.
* `VIBRANT_log_<input_file>`: this folder contains the log summary and run information.
