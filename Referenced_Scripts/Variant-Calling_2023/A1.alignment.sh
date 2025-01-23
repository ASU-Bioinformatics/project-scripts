#!/bin/bash

##### bwa alignment of sample strains to reference genome #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err               # STDERR (%j = JobId)
#SBATCH -t 0-12:00                    # estimated time needed
#SBATCH --mem=32G

module purge
module load mamba/latest

use="DIRECTORY"
help="FALSE"

VALID_ARGS=$(getopt -o a:r:g:f:l:h \
                    --long alignmentDir:,refDir:,refID:,fastqDir:,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --alignmentDir)
        echo "Alignment files will be written to the directory '$2'"
        alignmentDir="$2"
        shift 2
        ;;
    -r | --refDir)
        echo "The reference genome should be present in the directory '$2'"
        refDir="$2"
        shift 2
        ;;
    -g | --refID)
        echo "The reference genome and indexes have the prefix '$2'"
        refID="$2"
        shift 2
        ;;
    -f | --fastqDir)
        echo "FastQ files are in the directory '$2'"
        fastqDir="$2"
        shift 2
        ;;
    -l | --list)
        echo "Fastq sample IDs from '$2' will be analyzed"
        use="LIST"
        list="$2"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
        shift;
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
  This script aligns short-read sequences from a single species to the reference genome for that species.

  usage: sbatch A1.alignment.sh \
                -a /path/to/alignment/directory \
                -r /path/to/reference/directory \
                -g referenceID \
                -f /path/to/input/fastq/directory \
                (-l "sid1 sid2 sid3 ...") (-h)

  options:
    [ -a  |   --alignmentDir      |   directory to write alignment output files                                   ]
    [ -r  |   --refDir            |   directory containing the reference genome                                   ]
    [ -g  |   --refID             |   prefix for reference genome and indexes (will be created if not present)    ]
    [ -f  |   --fastqDir          |   directory containing the input fastq files for variant calling              ]
    [ -l  |   --list              |   optional; list of sample ids within fastqDir to use for alignment workflow  ]
    [ -h  |   --help              |   prints an informational message and exits script                            ]
EOF
  exit;
fi

source activate /data/biocore/programs/conda-envs/align_breakdance

if ! [[ -f "$refDir"/"$refID".fai ]]; then
  echo "Creating samtools faidx index"
  samtools faidx "$refDir"/"$refID"
fi

if ! [[ -f "$refDir"/"$refID".dict ]]; then
  echo "Creating picard sequence dictionary"
  picard CreateSequenceDictionary R="$refDir"/"$refID" O="$refDir"/"${refID%.*}".dict
fi

module load bwa-0.7.17-gcc-12.1.0

if ! [[ -f "$refDir"/"$refID".bwt ]]; then
  echo "Creating bwa index"
  bwa index "$refDir"/"$refID"
fi

mkdir -p "$alignmentDir"

cd "$alignmentDir"

if [[ $use == "DIRECTORY" ]]
then
  for i in $(find "$fastqDir" -maxdepth 1 -type f -name "*.fastq.gz" \
  | while read F; do basename $F; done \
  | cut -d "_" -f 1 | sort | uniq)
  do
    echo "Aligning $i"
    bwa mem -M -t 4 "$refDir"/"$refID" "$fastqDir"/"$i"_*_R1_001.fastq.gz "$fastqDir"/"$i"_*_R2_001.fastq.gz > "$i".sam
    echo "Converting $i to bam"
    samtools view -bS -o "$i".bam "$i".sam
    echo "Sorting $i bam"
    samtools sort -o "$i".sorted.bam "$i".bam
    echo "$i read groups"
    picard AddOrReplaceReadGroups I="$i".sorted.bam O="$i".sorted.rdgrp.bam \
                                  RGID=1 RGLB=ctrl RGPL=illumina RGPU=unit1 RGSM="$i".sorted \
                                  SO=coordinate VALIDATION_STRINGENCY=SILENT
    echo "indexing $i bams"
    picard BuildBamIndex INPUT="$i".sorted.rdgrp.bam OUTPUT="$i".sorted.rdgrp.bam.bai
    echo "checking alignment stats for $i"
    samtools flagstat "$i".sorted.rdgrp.bam > "$i".flagstat.txt
  done;

elif [[ $use == "LIST" ]]
then
  for i in "$list"
  do
    echo "Aligning $i"
    bwa mem -M -t 4 "$refDir"/"$refID" "$fastqDir"/"$i"_*_R1_001.fastq.gz "$fastqDir"/"$i"_*_R2_001.fastq.gz > "$i".sam
    echo "Converting $i to bam"
    samtools view -bS -o "$i".bam "$i".sam
    echo "Sorting $i bam"
    samtools sort -o "$i".sorted.bam "$i".bam
    echo "$i read groups"
    picard AddOrReplaceReadGroups I="$i".sorted.bam O="$i".sorted.rdgrp.bam \
                                  RGID=1 RGLB=ctrl RGPL=illumina RGPU=unit1 RGSM="$i".sorted \
                                  SO=coordinate VALIDATION_STRINGENCY=SILENT
    echo "indexing $i bams"
    picard BuildBamIndex INPUT="$i".sorted.rdgrp.bam OUTPUT="$i".sorted.rdgrp.bam.bai
    echo "checking alignment stats for $i"
    samtools flagstat "$i".sorted.rdgrp.bam > "$i".flagstat.txt
done;

fi;

#check flagstats for each sample before continuing to next step
