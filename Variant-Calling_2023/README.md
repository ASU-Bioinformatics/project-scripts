# Variant Calling Workflow

The goal of this pipeline is to identify SNPs and structural variants within one or more genomes as compared to a reference genome. After variants are detected, the impact of each variant is predicted and gene-specific annotations are generated for each one.

Additional Components to Develop:
* Functional annotations for genes with predicted High impact, for evaluation of functional enrichment.
* Genotype variant calling as an alternative to GATK's HaplotypeCaller (probably a module C)

## Module A: Alignment to Reference

The tools needed for this script are:
* [BWA]("https://doi.org/10.48550/arXiv.1303.3997", "Title")
* [SAMtools]("10.1093/gigascience/giab008", "Title")
* [picard]("http://broadinstitute.github.io/picard/", "Title")

These tools are installed in the biocore conda environment /data/biocore/programs/conda-envs/align_breakdance, which is hardcoded into the script. However, they can be accessed through Sol with ``module load bwa-0.7.17-gcc-12.1.0`` for BWA, ``module load picard-2.26.2-gcc-12.1.0`` for picard, and ``module load samtools-1.16-gcc-11.2.0`` for SAMtools. The activation line for the conda environment can be replaced with these three module loading commands to make the script usable for individuals without biocore access.

The primary inputs for this module are the reference genome and your short-read sequences of interest. [BWA]("https://doi.org/10.48550/arXiv.1303.3997", "Title") is used to align the short-read sequences to the reference genome, creating the FAI index file and the DICT dictionary file required by BWA if they aren't already present in the reference directory with the same prefix as the reference genome. In addition to the SAM/BAM output files, a flagstats quality metric files is created for each sample and should be checked before continuing with variant calling.

#### Example usage:

    sbatch A1.alignment.sh \
      -a /path/to/alignment/directory \
      -r /path/to/reference/directory \
      -g referenceID \
      -f /path/to/input/fastq/directory \
      (-l "sid1 sid2 sid3 ...") (-h)

#### Input Parameters:
  * ``-a,--alignmentDir``: directory to write alignment output files
  * ``-r,--refDir``: directory containing the reference genome (in fasta format)
  * ``-g,--refID``: prefix for reference genome and indexes (FAI and DICT files) - the indexes will be created if not already present
  * ``-f,--fastqDir``: directory containing the input fastq files for variant calling; by default, all fastq files in the directory will be called.
  * ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow, instead of everything in the directory
  * ``-h,--help``: prints an informational message and exits script

## Module B: Variant Calling

### Script B1: Small Variant Calling (SNPs and Small InDels)

Small variant calling - variants less than 50bp long, typically - are detected using GATK's HaplotypeCaller program. Both complete unfiltered VCF files and quality-filtered VCF files are output by the script. Default filtering parameters for the VCF file are derived from this Evolution and Genomics [tutorial]("http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/", "Title") and *Identification of Substrain-Specific Mutations by Massively Parallel Whole-Genome Resequencing of Synechocystis sp. PCC6803* by [Kanesaki et al., 2011]("https://doi-org.ezproxy1.lib.asu.edu/10.1093/dnares/dsr042", "Title"). This script is appropriate for single sample comparison to a genome, or a multiple-sample comparison when genotyping information isn't needed.

#### Example usage:

    sbatch B1.gatk.sh \
      -i /path/to/input/alignment/files/ \
      -r /path/to/reference/directory/reference.fasta \
      -o /path/to/variant/calling/output \
      -p 'INFO/DP<10 || FS>60.0' \
      (-l "sid1 sid2 sid3 ...") (-h)

#### Input Parameters:
* ``-i,--inputDir``: directory containing the input alignment files
* ``-r,--ref``: reference genome in fasta format
* ``-o,--outDir``: directory to write output files
* ``-p,--params``: parameter string for filtering VCF; it must be single-quoted. The default is: ``'(FMT/AD[0:1])/(FMT/DP)<0.6 || INFO/DP<10 || FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-10.5 || QD<2.0 || ReadPosRankSum<-8.0'``
* ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow
* ``-h,--help``: prints an informational message and exits script

### Script B2: Mid-Range Variant Calling

GRIDSS2 focuses on slightly longer variants than HaplotypeCaller, in my experience mostly identifying 50-300bp long breakpoints, insertions, and deletions. In addition, it incorporates a default quality filter for its variant calls, so the script outputs both the complete VCF (including poor quality calls) and the filtered VCF (including only calls that pass the GRIDSS2 filter).

The tools needed for this script include:
* [BWA]("https://doi.org/10.48550/arXiv.1303.3997", "Title") - tested with version 0.7.17
* [SAMtools]("10.1093/gigascience/giab008", "Title") - tested with version 1.16
* [OpenJDK]("https://openjdk.org/projects/jdk/12/", "Title") - version 11 or higher
* [htslib]("https://openjdk.org/projects/jdk/12/", "Title") - tested with version 1.16
* [R]("https://www.r-project.org/", "Title") - tested with version 4.2.2
* [VCFtools]("https://vcftools.github.io/index.html", "Title") - tested with version 0.1.14
* [GRIDSS2]("https://github.com/PapenfussLab/gridss", "Title") - tested with version 2.13.2

All but GRIDSS2 can be accessed through pre-installed Sol modules for ASU researchers, with the following commands:

    module load r-4.2.2-gcc-11.2.0
    module load htslib-1.16-gcc-11.2.0
    module load samtools-1.16-gcc-11.2.0
    module load jdk-12.0.2_10-gcc-12.1.0
    module load bwa-0.7.17-gcc-12.1.0
    module load vcftools-0.1.14-gcc-11.2.0

#### Example usage:

    sbatch B2.gridss2.sh \
     -i /path/to/input/alignment/files/ \
     -r /path/to/reference/directory/reference.fasta \
     -o /path/to/variant/calling/output \
     -s /path/to/directory/with/GRIDSS2/executables/ \
     (-l "sid1 sid2 sid3 ...") (-h)

#### Input Parameters:
* ``-i,--inputDir``: directory containing the input alignment files
* ``-r,--ref``: reference genome in fasta format, ideally with pathname
* ``-o,--outDir``: directory to write output VCF files
* ``-s,--scriptDir``: directory containing the GRIDSS2 executables
* ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow
* ``-h,--help``: prints an informational message and exits script

### Script B3: Long Structural Variant Identification

This script uses the program [DELLY]("https://github.com/dellytools/delly?tab=readme-ov-file", "Title") to identify long structural variants such as breakpoints and major insertions/deletions. Like GRIDSS2, DELLY incorporates its own quality filter so both unfiltered and filtered VCF files are output. The default output format is actually BCF, but this script includes a file type conversion to match the output from the other two variant calling scripts.

The tools needed for this script include:
* [VCFtools]("https://vcftools.github.io/index.html", "Title") - tested with version 0.1.14
* [BCFtools]("https://samtools.github.io/bcftools/", "Title") - tested with version 1.10.2
* [DELLY]("https://github.com/dellytools/delly?tab=readme-ov-file", "Title") - tested with version 1.2.6

VCFtools and BCFtools can both be accessed through Sol modules via ``module load vcftools-0.1.14-gcc-11.2.0`` and ``module load bcftools-1.10.2-gcc-12.1.0``. DELLY needs to be installed locally, ideally in a conda environment (which allows for a very straightforward installation and takes care of any Python dependencies).

#### Example usage:

    sbatch B3.delly.sh \
     -i /path/to/input/alignment/files/ \
     -r /path/to/reference/directory/reference.fasta \
     -o /path/to/variant/calling/output \
     (-l "sid1 sid2 sid3 ...") (-h)

#### Input Parameters:
* ``-i,--inputDir``: directory containing the input alignment files
* ``-r,--ref``: reference genome in fasta format, ideally with pathname
* ``-o,--outDir``: directory to write output BCF/VCF files
* ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow
* ``-h,--help``: prints an informational message and exits script

## Module C: snpEff Annotation

### Script C1: Standard snpEff Impact Prediction

This script uses [snpEff]("https://pcingola.github.io/SnpEff/", "Title") to predict the impact each variant may have on phenotype, based on the type of variant and genomic region in which it occurs (for example, a stop codon in the middle of a coding region will probably have a higher impact than a point mutation in an intron). It uses pre-made databases for the species of interest to identify the types of features impacted by each variant. snpEff will need to be locally installed (for biocore users, the core jar file is located in ``/data/biocore/programs/snpEff``). It also requires Java 11 or higher; on Sol, you can load OpenJDK for this requirement with ``module load jdk-12.0.2_10-gcc-12.1.0``.

#### Example usage:

    sbatch C1.snpeff-annotation.sh \
      -i /path/to/input/vcf/files/ \
      -r reference_species \
      -o /path/to/variant/calling/output/ \
      -s /path/to/snpEff/
      (-l "sid1 sid2 sid3 ...") (-h)

#### Input Parameters:
* ``-i,--inputDir``: directory containing the input VCF files
* ``-r,--ref``: reference species, using the name of the snpEff database in the snpEff/data directory
* ``-o,--outDir``: directory to write output annotated VCF files
* ``-s,--scriptDir``: directory containing the snpEff.jar file as well as the ./data directory containing the snpEff databases
* ``-l,--list``: optional; list of sample ids within fastqDir to use for alignment workflow
* ``-h,--help``: prints an informational message and exits script
