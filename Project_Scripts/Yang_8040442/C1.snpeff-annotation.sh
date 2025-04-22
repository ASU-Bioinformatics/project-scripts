#!/bin/bash

##### snpEff #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.snpEff.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.snpEff.err               # STDERR (%j = JobId)
#SBATCH -t 0-02:00                    # estimated time needed
#SBATCH -c 1
#SBATCH --mem=64G

umask 0007
module load jdk-12.0.2_10-gcc-12.1.0

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
        echo "All VCF files found in the directory '$2' will be analyzed."
        inputDir="$2"
        shift 2
        ;;
    -r | --ref)
        echo "The reference genome used is '$2'"
        ref="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Annotated VCF results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -s | --scriptDir)
        echo "The snpEff.jar file is found in the directory '$2'"
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
  This script runs snpEff to predict functional impact of variants.
  As input, it requires a VCF file and the name of the snpEff database stored in the snpEff directory

  usage: sbatch C1.snpeff-annotation.sh \
                  -i /path/to/input/vcf/files/ \
                  -r reference_species \
                  -o /path/to/variant/calling/output \
                  -s /path/to/snpEff \
                  (-l "sid1 sid2 sid3 ...") (-h)

  options:
    [ -i  |   --inputDir  |  directory containing the input VCF files                                    ]
    [ -r  |   --ref       |  snpEff database name for reference species                                  ]                                                                    ]
    [ -o  |   --outDir    |  directory to write output files                                             ]
    [ -s  |   --scriptDir |  directory containing the snpEff.jar file                                    ]
    [ -l  |   --list      |  optional; list of sample ids within fastqDir to use for alignment workflow  ]
    [ -h  |   --help      |  prints an informational message and exits script                            ]
EOF
  exit;
fi

if [ "$use" == "DIRECTORY" ];
then
  list=$(find "$inputDir" -type f -name "*.vcf" | \
          while read F; do basename $F; done | sort | uniq)
fi

cd $inputDir

mkdir -p "$outDir"

for i in $list
do
  java -Xmx16g -jar "$scriptDir"/snpEff.jar \
    "$ref" "$inputDir"/"$i" \
    -s "$outDir"/"${i%.vcf}".snpEff.summary.html > "$outDir"/"${i%.vcf}".snpeff.vcf
done
