#!/bin/bash

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%A.alignwrap.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.alignwrap.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-0:15
#SBATCH -c 1
#SBATCH --mem=2G

module purge

paired="TRUE"
use="DIRECTORY"
help="FALSE"

VALID_ARGS=$(getopt -o f:i:p:r:a:q:s:tl:h \
                    --long fastqDir:,intermediateDir:,starParams:,refDir:,alignmentOutDir:,quantOutDir:,scriptDir:,readTypeSingle,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --fastqDir)
        echo "The samples for alignment can be found in '$2'"
        fastqDir="$2"
        shift 2
        ;;
    -i | --intermediateDir)
        echo "FIFO server drive for intermediate alignment is '$2'"
        intermediateDir="$2"
        shift 2
        ;;
    -p | --starParams)
        echo "STAR alignment parameter text file is '$2'"
        starParams="$2"
        shift 2
        ;;
    -r | --refDir)
        echo "Reference files and STAR indexes are in '$2'"
        refDir="$2"
        shift 2
        ;;
    -a | --alignmentOutDir)
        echo "STAR alignment results will be output to the directory '$2'"
        alignmentOutDir="$2"
        shift 2
        ;;
    -q | --quantOutDir)
        echo "Stringtie alignment results will be output to the directory '$2'"
        quantOutDir="$2"
        shift 2
        ;;
    -s | --scriptDir)
        echo "Auxiliary scripts for alignment can be found in '$2'"
        scriptDir="$2"
        shift 2
        ;;
    -t | --readTypeSingle)
        echo "Samples are from single-end read sequencing"
        paired="FALSE"
        shift
        ;;
    -l | --list)
        echo "The list of samples in '$2' will be used instead of the complete fastq directory"
        use="LIST"
        list="$2"
        shift
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
  This script runs the whole alignment workflow, from STAR alignment through to the gene count matrix.
  All scripts A1 through A7 are included, and should be located in the directory provided by --scriptDir.

  usage: bash 2.alignment-wrapper.sh \
                -f /path/to/fastqs \
                -i /path/to/FIFO/intermediate \
                -p /path/to/starParams.txt \
                -r /path/to/reference/directory \
                -a /path/to/alignment/output \
                -q /path/to/stringtie/output \
                -s /path/to/auxiliary/scripts \
                (-l "sid1 sid2 sid3 ...")
                (-t) (-d) (-h)

  options:
    [ -f  |   --fastqDir          |   directory containing gzipped fastq files, named with standard Illumina formatting (ie, sid-1_S01_L001_R1.fastq*)  ]
    [ -i  |   --intermediateDir   |   FIFO-enabled server for STAR alignment                                                                              ]
    [ -p  |   --starParams        |   text file containing STAR parameters (see examples in scripts folder)                                               ]
    [ -r  |   --refDir            |   directory containing the reference genome, GTF file, and STAR indexes - must contain STAR genomeParamers.txt        ]
    [ -a  |   --alignmentOutDir   |   location for output alignment files (eg, sorted bam files)                                                          ]
    [ -q  |   --quantOutDir       |   location for output stringtie quantification files (eg, transcript abundance tables for each sample)                ]
    [ -s  |   --scriptDir         |   directory where the auxiliary python and perl scripts can be found (scripts A4-A7)                                  ]
    [ -t  |   --readTypeSingle    |   optional; specifies that input reads are single-end (unidirectional)                                                ]
    [ -l  |   --list              |   optional; list of sample ids within fastqDir to use for alignment workflow                                          ]
    [ -h  |   --help              |   prints an informational message and exits script                                                                    ]
EOF
  exit;
fi

pathGTF=$(grep 'sjdbGTFfile' "$refDir"/genomeParameters.txt | tail -n 1 | awk '{ printf $2 }')
refGTF=$(basename $pathGTF)

if [ "$use" == "DIRECTORY" ];
then
  list=$(find "$fastqDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
fi

#create array of sample IDs and count sample number
samples=($list)
count=${#samples[@]}
echo $count
echo $list

align=$(sbatch --array=0-$((count-1))%6 --parsable --export=ALL \
        "$scriptDir"/A1.cpu-align.sh \
          --array "$list" \
          --fastqDir "$fastqDir" \
          --intermediateDir "$intermediateDir" \
          --starParams "$starParams" \
          --refDir "$refDir" \
          --outputDir "$alignmentOutDir")

echo $align

quant=$(sbatch --dependency=afterok:"$align" --array=0-$((count-1)) --parsable --export=ALL \
          "$scriptDir"/A2.stringtie-quant.sh \
            --list "$list" \
            --alignmentDir "$alignmentOutDir" \
            --quantDir "$quantOutDir" \
            --refGTF "$refDir"/"$refGTF")

merge=$(sbatch --dependency=afterok:"$quant" --parsable --export=ALL \
        "$scriptDir"/A3.merge-quant.sh \
          --scriptDir "$scriptDir" \
          --quantDir "$quantOutDir" \
          --alignmentDir "$alignmentOutDir")
