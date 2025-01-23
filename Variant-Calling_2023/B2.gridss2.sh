#!/bin/bash

##### large structural variant detection #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.gridss2.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.gridss2.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00                    # estimated time needed
#SBATCH -c 8
#SBATCH --mem=40G

module purge
module load r-4.2.2-gcc-11.2.0
module load htslib-1.16-gcc-11.2.0
module load samtools-1.16-gcc-11.2.0
module load jdk-12.0.2_10-gcc-12.1.0
module load bwa-0.7.17-gcc-12.1.0
module load vcftools-0.1.14-gcc-11.2.0

use="DIRECTORY"
help="FALSE"

VALID_ARGS=$(getopt -o i:r:o:s:l:h \
                    --long help,inputDir:,ref:,outDir:,scriptDir:,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "All alignment files found in the directory '$2' will be analyzed."
        inputDir="$2"
        shift 2
        ;;
    -r | --ref)
        echo "The reference genome used is '$2'"
        ref="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "VCF results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -s | --scriptDir)
        echo "GRIDSS2 executables are found in the directory '$2'"
        scriptDir="$2"
        shift 2
        ;;
    -l | --list)
        echo "Alignment files with sample IDs in '$2' will be analyzed"
        use="LIST"
        list="$2"
        shift 2
        ;;
    -h | --help )
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
  This script runs GRIDSS2 to identify medium length (32-400bp) structural variants in WGS samples.
  As input, it requires a sorted BAM file with read groups added and duplicates marked for each sample
  * The suffix for these BAM files should be '.sorted.rdgrp.mrkd.bam'.
  * Everything before the first '.' in the sample name will be used as the sample ID
  Additionally, it requires the reference genome fasta file to which they were aligned
  All samples from the alignment folder will be used unless the -l option provides a quoted, space-separated list of sample ids
  GRIDSS2 incorporates its own quality filter, so an unfiltered and filtered VCF files are output.
  This script is NOT suitable for genotype generation/comparison of multiple samples.

  usage: sbatch B2.gridss2.sh \
          -i /path/to/input/alignment/files/ \
          -r /path/to/reference/directory/reference.fasta \
          -o /path/to/variant/calling/output \
          -s /path/to/directory/with/GRIDSS2/executables/ \
          (-l "sid1 sid2 sid3 ...") (-h)

  options:
    [ -i  |   --inputDir     |  directory containing the input alignment files                              ]
    [ -r  |   --ref       |  the reference genome fasta file, ideally with pathname                      ]                                                                    ]
    [ -o  |   --outDir    |  directory to write output files                                             ]
    [ -s  |   --scriptDir |  directory containing the GRIDSS2 executables                                ]
    [ -l  |   --list      |  optional; list of sample ids within fastqDir to use for alignment workflow  ]
    [ -h  |   --help      |  prints an informational message and exits script                            ]
EOF
  exit;
fi

if [ "$use" == "DIRECTORY" ];
then
  list=$(find "$inputDir" -type f -name "*.sorted.rdgrp.mrkd.bam" | \
          while read F; do basename $F; done | cut -d '.' -f 1 | sort | uniq)
fi

for i in $list
do
  $scriptDir/scripts/gridss \
    -t 8 -r $ref \
    -j $scriptDir/gridss-2.13.2-gridss-jar-with-dependencies.jar \
    -o $outDir/"$i".gridss.vcf \
    $inputDir/"$i".sorted.rdgrp.mrkd.bam

  vcftools --vcf $outDir/"$i".gridss.vcf \
           --out $outDir/"$i".gridss.hardfiltered.vcf \
           --remove-filtered-all --recode
done
