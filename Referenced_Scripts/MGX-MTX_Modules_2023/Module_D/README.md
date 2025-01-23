# Module D
This module creates a contig assembly integrating all reads from the dataset; this strategy allows for more complete MAGs and more thorough alignment and abundance profiling, especially from species that are underrepresented in a specific sample but more prevalent in another sample that can then provide a reference for the downstream alignment.

This module is only for MGX and unpaired MTX data. Paired MTX data uses the assembly from the corresponding MGX data for transcript alignment and profiling.

On Sol, the scripts in this module use the environments `/data/biocore/programs/conda-envs/megahit-args`, `/data/biocore/programs/conda-envs/quast-env`, ``/data/biocore/programs/conda-envs/bowtie2-env``, ``/data/biocore/programs/mamba-envs/gfastats-env``, and ``/data/biocore/programs/conda-envs/bb-env``. (Note to self, try to streamline these environments!) The environments are hardcoded into the scripts, but if you are not able to access the installations, you can use the included yml files to create equivalent environments and replace the environment pathname in the script.

## D1.megahit-assemble.sh
This script uses the [MEGAHIT](https://github.com/voutcn/megahit,"Title") assembler to create a contig file, then runs [QUAST](https://github.com/ablab/quast,"Title") to obtain basic stats about the contigs in the assembly.

On Sol, this script uses the environment ``/data/biocore/programs/conda-envs/megahit-args`` for running MEGAHIT and the environment ``/data/biocore/programs/conda-envs/quast-env`` for running QUAST.

A sample call:

    sbatch ./D1.megahit-assemble.sh \
      --inputDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly

The arguments that can be passed are as follows:
* ``-i,--inputDir``: the folder containing the preprocessed fastq.gz files
* ``-o,--outDir``: the folder to which MEGAHIT should write output files. This directory cannot exist before the script is one or MEGAHIT will throw an error, unless the ``--resume`` option is used to finish building an incomplete assembly.
* ``-r,--resume``: this flag should be included to finish building a partial assembly. With the 2 day limit on Sol for high memory nodes, it almost always times out before completing the assembly and will need a second call to complete.

Output files from MEGAHIT include the assembly itself, named ``final.contigs.fa`` as well as a folder containing all the intermediate files produced by the workflow. From the QUAST run there will be a folder containing general statistical information about the contigs (length, number, and so on).

## D2.bowtie-index.sh
This script builds a [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml,"Title") index for the contig assembly (hardcoded as the MEGAHIT default ``final.contigs.fa``) to allow for sample alignment back to the assembly in the next step. It also uses [gfastats](https://github.com/vgl-hub/gfastats,"Title") to create an assembly graph, although the statistic reports don't appear to be working correctly.

On Sol, Conda environments for Bowtie2 and gfastats are available at ``/data/biocore/programs/conda-envs/bowtie2-env`` and ``/data/biocore/programs/mamba-envs/gfastats-env`` respectively, and are hardcoded into the script.

A sample call:

    sbatch ./D2.bowtie-index.sh \
      --fastqDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --assemblyDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly \
      --prefix coassembly

The script takes the following arguments:
* ``-f,--fastqDir``: the directory containing the processed fastq.gz files
* ``-a,--assemblyDir``: the directory containing the contig coassembly from D1
* ``-p,--prefix``: the prefix for the bowtie2 index files

Output files are bowtie2 index files (either small or large, automatically determined by Bowtie2) for the assembly as well as a GFA assembly graph file and corresponding statistic report text files.

## D3.bowtie-align.wrapper.sh
This script runs Bowtie2 alignment for all samples in a list or directory to the assembly fasta file from D1 and the Bowtie2 index from D2. Alignment characterizes the relative abundance of each contig in the assembly in the reads for each sample and assists with generating MAG bins in module G (as the binning tool can assume that all pieces of the same genome have approximately the same read depth in any given sample).

In addition, the script converts the SAM files to sorted BAM file and calculates coverage statistics and RPKM values for each alignment, using [samtools](https://www.htslib.org/,"Title") and [BBMap](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/bb-tools-user-guide/bbmap-guide/,"Title")'s pileup tool.

As with D2, the Conda environment ``/data/biocore/programs/conda-envs/bowtie2-env/`` is hardcoded into the script for use on Sol. If this environment variable is changed, please also change the PERL5LIB environment variable to specifically use the installation of Perl within the Conda environment. Samtools appears to be natively available on Sol, and BBMap can be accessed through the Conda environment ``/data/biocore/programs/conda-envs/bb-env``.

A sample call:

    bash ./D3.bowtie-align.wrapper.sh \
      --fastqDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/fastq \
      --assembly /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/megahit-assembly/coassembly \
      --outDir /data/gencore/analysis_projects/6209010_Keaton_Meta/DNA-Metagenomics/human-removed/bowtie2-alignments

The arguments for this script are:
* `-f,--fastqDir`: the full pathname to the directory where the fastq.gz files are stored
* ``-a,--assembly``: the full pathname to the bowtie2 index for the assembly, including the index prefix
* ``-o,--outDir``: the full pathname to the directory where out files should be written
* ``-l,--list``: a list of sample IDs in the fastq directory that should be aligned; this is optional and not required if all samples in the directory need to be aligned.

Output files include SAM files, SAM files with a problematic header section removed, BAM files, sorted BAM files, and BAI index files for the sorted BAM files. Additionally, a covstats and rpkm text file are created for each sample.

## Before moving on
Check the QUAST output for the assembly to see if the contig length makes sense (for example, MTX contigs should be fairly short overall).
