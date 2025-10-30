#!/bin/bash

##### qiime2 pt1: dada2 #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 3-0:00                         # estimated time needed (dada2 can take a while)
#SBATCH -c 2
#SBATCH --mem=128G

umask 0007

fastqDir="$pwd"/fastq
qiimeDir="$pwd"/qiime2
metadata="$pwd"/metadata.txt
environment="/data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7/"
inputPairing="paired"
dada2="single"
cutLen=0
manifest="FALSE"
help="FALSE"
runDemux="FALSE"
runDada2="FALSE"
runStats="FALSE"
mode="all"

VALID_ARGS=$(getopt -o f:q:m:p:i:o:e:h \
                    --long fastqDir:,qiimeDir:,metadata:,pairing:,inputManifest:,mode:,environment:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --fastqDir)
        echo "Input fastq files will be read from '$2'"
        fastqDir="$2"
        shift 2
        ;;
    -q | --qiimeDir)
        echo "Output Qiime2 files will be written to '$2'"
        qiimeDir="$2"
        shift 2
        ;;
    -m | --metadata)
        echo "The sample metadata file to use is '$2'"
        metadata="$2"
        shift 2
        ;;
    -p | --pairing)
        if [ "$2" == "p" ];
          then
            echo "The read type is paired-end"
            inputPairing="paired"
            dada2="paired"
          elif [ "$2" == "s" ];
          then
            echo "The read type is single-end"
            inputPairing="single"
            dada2="single"
          elif [ "$2" == "ps" ];
          then
            echo "The read type is paired-end but DADA2 will be run with single-end parameters"
            inputPairing="paired"
            dada2="single"
        fi

        shift 2
        ;;
    -i | --inputManifest)
        echo "The absolute filepaths of the fastq files to evaluate can be found in '$2'"
        manifest="$2"
        shift 2
        ;;
    -o | --mode)
        mode="$2"
        for x in $mode;
        do
          if [ $x == "all" ];
          then
            runDemux="TRUE"
            runDada2="TRUE"
            runStats="TRUE"
          fi
          if [ $x == "demux" ];
          then
            runDemux="TRUE"
          fi
          if [ $x == "dada2" ];
          then
            runDada2="TRUE"
          fi
          if [ $x == "stats" ];
          then
            runStats="TRUE"
          fi
        done
        echo "The following steps will be performed: '$2'"
        shift 2
        ;;
    -e | --environment)
        echo "The conda environment to use is '$2'"
        environment="$2"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
        break
        ;;
    --)
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1 - please correct."
        ;;
  esac
done

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script runs the first part of our standard Qiime2 analysis (loading the data, denoising,
  and summarizing statistics by metadata information).

  usage: sbatch qiime2_pt1_args.sh
            -f /path/to/fastq-input -q /path/to/qiime-output
            -m /path/to/metadata.txt -p pairing
             -e /path/to/conda/environment (-h)

  options:
    [ -f  |   --fastqDir     |   directory containing gzipped fastq files, named with standard Illumina formatting (ie, sid-1_S01_L001_R1.fastq*)               ]
    [ -q  |   --qiimeDir     |   directory for Qiime2 output files (will be created if it doesn't already exist; previous files will be overwritten)            ]
    [ -m  |   --metadata     |   text file containing metadata information in Qiime2-compatible format                                                          ]
    [ -p  |   --pairing      |   the read pairing for the sequencing. Allowable options are p, s, and ps (paired, single, and paired load with single DADA2)    ]
    [ -i  |   --inputFiles   |   the absolute path of a fastq manifest file in Qiime2-compatible format; default is to generate this from the fastqDir          ]
    [ -o  |   --mode         |   the modules within the script to run - options are "all" or any combination of "demux", "dada2", and "stats" in a quoted list  ]
                                   to run "dada2", input demux files must be in the Qiime2 output folder and named according to the script parameters for       ]
                                   paired-end or single-end reads (demux-paired.qza or demux-single.qza)                                                        ]
                                   to run "stats", input dada2 files must be in the Qiime2 output folder and named according to the script parameters for       ]
                                   read pairing and dada2 pairing; for example, stats-dada2-paired-single.qza for a file generated with paired-end input reads  ]
                                   and single-end DADA2 processing. The stats, table, and rep-seqs qza files need to be present in the Qiime2 output directory  ]
    [ -e  |   --environment  |   location for the Qiime2 environment to activate                                                                                ]
    [ -h  |   --help         |   prints an informational message and exits script                                                                               ]
EOF
  exit;
fi

module purge

mkdir -p "$qiimeDir"

module load mamba/latest
source activate "$environment"

# make fastq manifest file - this lets us run qiime without having to rename the files
# if they aren't in the traditional Casava naming format (which the more modern machines aren't)

if [ "$manifest" == "FALSE" ];
then
  cd $fastqDir

  touch int1.txt
  echo "sample-id" > int1.txt
  find $(readlink -f ..) -type f -name "*_R2*.fastq*" | sort | uniq \
    | while read F; do basename $F; done \
    |cut -d '_' -f 1 >> int1.txt

  touch int2.txt

  if [ "$pairing" == "s" ];
  then
    echo "absolute-filepath" > int2.txt
    for i in $(find $(readlink -f ..) -type f -name "*_R1*.fastq*" | sort | uniq ); do echo $i; done >> int2.txt
    paste int1.txt int2.txt > fq_manifest.tsv
    rm int*.txt
  else
    echo "forward-absolute-filepath" > int2.txt
    for i in $(find $(readlink -f ..) -type f -name "*_R1*.fastq*" | sort | uniq ); do echo $i; done >> int2.txt
    touch int3.txt

    echo "reverse-absolute-filepath" > int3.txt
    for i in $(find $(readlink -f ..) -type f -name "*_R2*.fastq*" | sort | uniq ); do echo $i; done >> int3.txt
    paste int1.txt int2.txt int3.txt > fq_manifest.tsv
    rm int*.txt
  fi

  manifest="$fastqDir"/fq_manifest.tsv
fi

cd "$qiimeDir"

# if specified mode is "all" or contains "demux", the script will import the data using the fastq manifest document
# next, it will create a visual summary of the demultiplexed data
if [ "$runDemux" == "TRUE" ];
then
  if [ "$pairing" == "s" ];
  then
    qiime tools import \
    --type 'SampleData[SequencesWithQuality]' \
    --input-path "$manifest" \
    --input-format SingleEndFastqManifestPhred33V2 \
    --output-path "$qiimeDir"/demux-"$inputPairing".qza
  else
    qiime tools import \
      --type 'SampleData[PairedEndSequencesWithQuality]' \
      --input-path "$manifest" \
      --input-format PairedEndFastqManifestPhred33V2 \
      --output-path "$qiimeDir"/demux-"$inputPairing".qza
  fi

  qiime demux summarize \
    --i-data "$qiimeDir"/demux-"$inputPairing".qza \
    --o-visualization "$qiimeDir"/demux-"$inputPairing".qzv
fi

# if specified mode is "all" or contains "dada2", the script will use the demux files found in the qiime2 output folder to run data2 filtering
# make sure that the input demux files are named according to the script parameters for paired-end or single-end reads (demux-paired.qza or demux-single.qza)
# dada2 will filter any phiX reads, low quality reads, and chimeric reads from the input data
if [ "$runDada2" == "TRUE" ];
then
  if [ "$dada2" == "single" ];
  then
    qiime dada2 denoise-single \
      --i-demultiplexed-seqs "$qiimeDir"/demux-"$inputPairing".qza \
      --p-trim-left 0 \
      --p-trunc-len $cutLen \
      --p-n-threads 0 \
      --o-representative-sequences "$qiimeDir"/rep-seqs-"$inputPairing"-"$dada2".qza \
      --o-table "$qiimeDir"/table-"$inputPairing"-"$dada2".qza \
      --o-denoising-stats "$qiimeDir"/stats-dada2-"$inputPairing"-"$dada2".qza
  else
    qiime dada2 denoise-paired \
      --i-demultiplexed-seqs "$qiimeDir"/demux-"$inputPairing".qza \
      --p-trunc-len-f $cutLen \
      --p-trunc-len-r $cutLen \
      --p-n-threads 0 \
      --o-representative-sequences "$qiimeDir"/rep-seqs-"$inputPairing"-"$dada2".qza \
      --o-table "$qiimeDir"/table-"$inputPairing"-"$dada2".qza \
      --o-denoising-stats "$qiimeDir"/stats-dada2-"$inputPairing"-"$dada2".qza
  fi
fi

# if the specified mode is "all" or contains "stats" run the metadata and feature table summarize and tabulate scripts
# this script requires the dada2 output Qiime2 artefacts, named according to both the input pairing and the DADA2 pairing
# for example, stats-dada2-paired-single.qza for a file generated with paired-end input reads and single-end DADA2 processing
# in addition to the stats qza file, the table qza file and rep-seqs qza file need to be present in the Qiime2 output directory
if [ "$runStats" == "TRUE" ];
then
  qiime metadata tabulate \
    --m-input-file "$qiimeDir"/stats-dada2-"$inputPairing"-"$dada2".qza \
    --o-visualization "$qiimeDir"/stats-dada2-"$inputPairing"-"$dada2".qzv

  qiime feature-table summarize \
    --i-table "$qiimeDir"/table-"$inputPairing"-"$dada2".qza \
    --o-visualization "$qiimeDir"/table-"$inputPairing"-"$dada2".qzv \
    --m-sample-metadata-file "$metadata"

  qiime feature-table tabulate-seqs \
    --i-data "$qiimeDir"/rep-seqs-"$inputPairing"-"$dada2".qza \
    --o-visualization "$qiimeDir"/rep-seqs-"$inputPairing"-"$dada2".qzv
fi

chmod -R g+w *
