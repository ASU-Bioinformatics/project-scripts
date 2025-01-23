# Module A

Module A contains two scripts (one of which is optional), but it's one of the most essential modules for the entire pathway because *de novo* assembly is such a large part of the analysis. The assembler will attempt to use any adapter sequences in the reads to assemble those reads into contigs, which is obviously going to lead to misassembled contigs!
Similarly, quality filtering helps reduce both misassembly and misalignment by removing sequencing errors that could muddy the differences between strains and species in the data.

## A1.trimmomatic.sh
 A1 requires only the directory containing the fastq.gz files as input. The trimmomatic call will by default remove all read pairs with an average Phred score of less than 15 in any 4bp window across both the forward and reverse read, and will trim adapter sequences off the ends of the reads using the trimmomatic TruSeq3-PE-2 adapter set.

A1 uses both [Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic,"Title") and [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/,"Title"), both of which can be loaded as built-in Sol modules. You can call the A1 script on Sol using the following command, without needing access to any specific conda environments:

     bash ./A1.trimmomatic.sh /data/gencore/analysis_projects/6078853_Otak_DNA/fastq/original

The A1 script ends after generating fastqc quality metric files for the filtered and trimmed reads. I typically will then run [MultiQC](https://multiqc.info/, "Title") manually to summarize these files. MultiQC requires Python <= 3.7 so it will most likely require a dedicated conda environment. For the bioinformatics core, this can be found at ``/data/biocore/programs/multiqc_py3.7/``; it can easily be duplicated on your system using the included yml file. On Sol, MultiQC needs to be run from your home directory for I/O reasons.

## A2.remove-host.sh
A2 only needs to be used for a host-based sample such as human fecal material; it is not necessary for environmental samples. If there is a risk of contamination it won't hurt, though.

This script requires [Bowtie2](https://github.com/BenLangmead/bowtie2,"Title"), [samtools](https://www.htslib.org/,"Title"), and pileup from the [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/,"Title") suite of tools (for alignment stats). On Sol, these are provided through conda environments on `/data/biocore`: `/data/biocore/programs/conda-envs/bowtie2-env/` and `/data/biocore/programs/conda-envs/bb-env`; they can be installed on your local environment using the provided yml files, after which the hardcoded environment pathnames in the scripts will need to be edited.

A sample call:

    bash ./A2.remove-host.wrapper.sh \
      --assembly /data/gencore/databases/reference_genomes/human/Homo_sapiens/Ensembl/GRCh37/Sequence/Bowtie2Index/genome \
      --fqDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/fastq/adapter-trimmed/fastq/ \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --use "DIRECTORY"


The options are as follows; only -l is optional:
* ``-a,--assembly``: full pathname to the reference genome bowtie2 index, including the index prefix
* ``-f,--fqDir``: full pathname to the input directory containing the fastq.gz files
* ``-o,--outDir``: full pathname to the output directory
* ``-u,--use``: use case (DIRECTORY or LIST)
* ``-l,--list``: list containing samples names, if LIST is chosen

The output will be fastq.gz files containing all read pairs that did not align concordantly to the host, SAM files containing alignment data for all reads, alignment statistic text files, and fastqc quality files. As with A1, I typically multiqc to summarize the quality and quantity of the unmapped reads.

## A3.remove-host-rna.sh
A3 is the RNA version of A2, and uses the [STAR aligner](https://github.com/alexdobin/STAR,"Title") instead of Bowtie2. STAR 2.7.10a is installed on `/data/biocore/` and can be run from the directory `/data/biocore/programs/STAR-2.7.10a/source`, although the script expects to find it in your path. Just make sure to use 2.7.10a STAR indexes for the reference.

STAR doesn't work well on the ``/data`` drives due to I/O issues, so you'll want to set a ``/scratch`` directory for the tempDir argument.

A sample call:

    bash ./A3.remove-host-rna.wrapper.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/fastq/adapter-trimmed/fastq \
      --refDir /data/gencore/databases/reference_genomes/human/Ensembl_GRCh38_102_agaveSTAR2.7.3/star_2.7.10a_indexes \
      --tempDir /scratch/kawoodbu \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/RNA-Metagenomics/human-removed \
      --use "DIRECTORY"

The options are as follows; only -l is optional:
* ``-r,--refdir``: full pathname to the directory containing the STAR indexes
* ``-i,--inputdir``: full pathname to the input directory containing the fastq.gz files
* ``-o,--outDir``: full pathname to the output directory
* ``-u,--use``: use case (DIRECTORY or LIST)
* ``-l,--l``: list containing samples names, if LIST is specified for ``--use`` option

## BEFORE MOVING ON

Check **every sample** to ensure that read depth and quality are sufficient. Liu et al. recommend at least 15 million reads in [their 2022 paper](https://doi.org/10.1139/gen-2021-0120,"Title"), so that's a good cutoff value when possible. Anything with less than 5 million reads should probably be removed from further analysis, and anything under 1 million reads definitely should be dropped.
