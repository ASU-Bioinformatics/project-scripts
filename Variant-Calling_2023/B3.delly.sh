#!/bin/bash

##### large structural variant detection #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.delly.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.delly.err               # STDERR (%j = JobId)
#SBATCH -t 0-08:00                    # estimated time needed
#SBATCH -c 1
#SBATCH --mem=64G

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/delly-env

module load vcftools-0.1.14-gcc-11.2.0
module load bcftools-1.10.2-gcc-12.1.0

use="DIRECTORY"
help="FALSE"

VALID_ARGS=$(getopt -o i:r:o:l:h \
                    --long inputDir:,ref:,outDir:,list:,help \
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
  This script runs DELLY to identify long (50-1kb) structural variants in WGS samples.
  As input, it requires a sorted BAM file with read groups added and duplicates marked for each sample
  * The suffix for these BAM files should be '.sorted.rdgrp.mrkd.bam'.
  * Everything before the first '.' in the sample name will be used as the sample ID
  Additionally, it requires the reference genome fasta file to which they were aligned
  All samples from the alignment folder will be used unless the -l option provides a quoted, space-separated list of sample ids
  DELLY incorporates its own quality filter, so unfiltered and filtered VCF files are output, along with the default BCF output.
  This script is NOT suitable for genotype generation/comparison of multiple samples.

  usage: sbatch B3.delly.sh \
          -i /path/to/input/alignment/files/ \
          -r /path/to/reference/directory/reference.fasta \
          -o /path/to/variant/calling/output \
          (-l "sid1 sid2 sid3 ...") (-h)

  options:
    [ -i  |   --inputDir   |  directory containing the input alignment files                                                                                ]
    [ -r  |   --ref       |  the reference genome fasta file, ideally with pathname                                                                                            ]
    [ -o  |   --outDir  |  directory to write output files                                                                                               ]
    [ -l  |   --list     |   optional; list of sample ids within fastqDir to use for alignment workflow                                                   ]
    [ -h  |   --help     |   prints an informational message and exits script                                                                             ]
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
  delly call -g $ref -o "$outDir"/"$i".delly.bcf "$inputDir"/"$i".sorted.rdgrp.mrkd.bam

  bcftools convert "$outDir"/"$i".delly.bcf -o "$outDir"/"$i".delly.vcf

  vcftools --vcf "$outDir"/"$i".delly.vcf \
           --out "$outDir"/"$i".delly.hardfiltered.vcf \
           --remove-filtered-all --recode
done
