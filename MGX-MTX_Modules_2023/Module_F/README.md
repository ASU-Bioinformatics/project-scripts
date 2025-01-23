# Module F
This module predicts proteins on the contigs from the coassembly, annotates them, and outputs classification data as well as an overview of the functional KEGG modules present in the coassembly. This module is only needed where there is a unique coassembly to use - that is, for MGX and standalone MTX data, but not paired MTX data since that will rely on the DNA coassembly file and its predicted proteins.

Different scripts in this module rely on pre-installed Conda environments on Sol. If you don't have access to the installation or are working on a different system, you can create the environments using the files `catbat.yml` and `microbeAnnotator.yml` and replace the environment activation paths within the code.

## F1.cat-contig-annotate.sh
This script uses [CAT/BAT](https://github.com/dutilh/CAT,"Title") to predict proteins on the coassembly contigs and determine taxonomic classification of the contigs using those predictions. It takes a CAT database as input (the prepared database 2021-01-07_CAT_database is available in the gencore data mount on Sol), along with the coassembly contig fasta file.

On Sol, this script requires access to the Conda environment `/data/biocore/programs/mamba-envs/catbat-env`, which is hardcoded into the script.

A sample call:

    sbatch ./F1.cat-contig-annotate.sh \
      --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
      --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
      --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

Input parameters:
* `-d,--databaseDir`: the directory containing the prepared CAT database. For the 2021-02-07 downloadable database, that includes:
    * 2021-01-07.nr.dmnd
    * 2021-01-07.nr.fastaid2LCAtaxid
    * 2021-01-07.nr.gz
    * 2021-01-07.nr.gz.md5
    * 2021-01-07.nr.taxids_with_multiple_offspring
* `-t,--taxaDir`: the directory containing the taxonomy files for the prepared CAT database. This should include the following files:
    * citations.dmp
    * division.dmp
    * gencode.dmp
    * names.dmp
    * delnodes.dmp
    * gc.prt
    * merged.dmp
    * nodes.dmp
* `-c,--contigs`: the coassembly contig file
* `-o,--outDir`: the directory into which the annotation files should be written

Output includes an amino acid fasta file containing the predicted protein sequences, a gff file mapping the predicted proteins onto the contigs, a summary file showing the number of hits and top lineage score for each contig, and a classification file detailing the taxID assigned to each contig (which is used to  make classification files containing human readable taxa names in script F4.)

## F2.microbe-annotator.sh and wrapper script
This script uses the tool [MicrobeAnnotator](https://github.com/cruizperez/MicrobeAnnotator,"Title") to functionally annotate the proteins predicted by CAT, primarily aiming to assign KO numbers to each protein by iteratively searching the databases [KOfamscan](https://github.com/takaram/kofam_scan,"Title"), [UniProt](https://www.uniprot.org/,"Title")/Swissprot, [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/,"Title"), and UniProt/Trembl.

On Sol, this script requires the conda environment `/data/biocore/programs/mamba-envs/microbe-annotator-env/`, which is hardcoded into the script.

A sample call:

    sbatch ./F2.microbe-annotator.wrapper.sh \
      --databaseDir /data/gencore/databases/microbe-annotator --split \
      --proteinFasta /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations/out.CAT.predicted_proteins.faa

Input parameters:
* `-d,--databaseDir`: the directory containing the MicrobeAnnotator database
* `-p,--proteinFasta`: the .faa file containing amino acid sequences for predicted proteins from the contig coassembly file
* `-s,--split`: calling this option splits the predicted protein file into 10 smaller files to run annotation on in parallel, to increase overall efficiency and avoid resource limits on Sol.
* `-c,--continue`: calling this option allows an incomplete run to be resumed. You do not use the `--split` option with this call if the split .faa files are still present from the original run!
* `-l,--list`: this option allows the user to use the `--continue` option with only a subset of the original split files, to avoid rerunning the annotation for the subset that did complete successfully.

Output files include a results folder for each reference database, a table, barplot, and heatmap for the metabolic summary data, and the folder `annotation_results` containing all functional annotations for each predicted protein and a list of all KO numbers identified in the assembly. The barplot and heatmap produced here are not particularly useful since there is only one sample (the coassembly). Scripts in Module H will create similar graphics for all the individual samples in comparison to one another.

## F3.microbe-annotator-merge.sh
This script is used to combine the annotations from .faa files split in F2. If you ran F2 without the `--split` option, you can skip this one. It only uses standard Unix commands, so no special tools or environments need to be installed.

A sample call:

    sbatch ./F3.microbe-annotator-merge.sh /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

Input parameters are positional and only include the directory that contains the split microbe-annotator output folders. This directory should contain the files in the locations they were created by scripts F1 and F2 as it relies on that file structure.

Output parameters include .annot and .ko files for the complete set of predicted proteins, both with (prefix 'final') and without (prefix 'all') a header row. This script will also create a new directory for unmerged annotations and move the files specific to the split protein files into it.

## F4.ko_mapper.sh
This script makes a few auxiliary files for CAT and MicrobeAnnotator that are useful for visualization or downstream analysis. First, it adds the taxa names to the CAT classification files (instead of just taxa IDs). Second, it runs MicrobeAnnotator's script ko_mapper.py to create a barplot, heatmap, and metabolic summary table for the merged data. While the figures here aren't very useful, the metabolic summary may be.

On Sol, this script uses the conda environments at `/data/biocore/programs/mamba-envs/catbat-env` and `/data/biocore/programs/mamba-envs/microbe-annotator-env`, which are hardcoded into the script. The taxonomy directory and script directory on Sol are also coded into the script as default values, but can be overridden with the `--taxDir` and `--scriptDir` options.

A sample call:

    sbatch ./F4.ko_mapper.sh \
      --taxDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
      --scriptDir /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline \
      --contigs /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/final.contigs.fa \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/contig-annotations

Input parameters:
* `-t,--taxDir`: this is the directory where the taxonomy files for the CAT database can be found
* `-s,--scriptDir`: this is the directory within the conda environment where MicrobeAnnotator is installed that contains the ko_mapper.py script and the other files referenced by that script. It won't work if you just copy the ko_mapper.py script on its own out to another folder!
* `-c,--contigs`: this is the contig assembly file that was used as input for the earlier scripts in this module
* `-o,--outDir`: this is the directory to which output files should be written, and is also the file in which the standard output from the earlier scripts in this module were written, as F4 relies on that default file structure.

## Before Moving On...
Take a look at the barplot and heatmap from the merged annotations. In my experience, there are so many complete functional modules that it isn't possible to even read the values! If you're only seeing a few complete modules, that would be a red flag about the success of your protein predictions in my opinion.
