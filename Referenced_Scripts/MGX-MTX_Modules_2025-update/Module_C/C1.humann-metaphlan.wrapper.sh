#!/bin/bash

##### humann3 test #####

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.out
#SBATCH -e slurm.%j.err
#SBATCH -t 0-12:00

#### notes on C1.humann-metaphlan ####
# this script takes a long time to run, will need to be run for each sample independently
# set $concatenate to TRUE if the fastq files need to be concatenated still
# set $concatenate to FALSE if the fastq files have already been concatenated
# set $type to DNA for metagenomics or metatranscriptomics without paired metagenomics data
# set $type to RNA for metatranscriptomics with a paired metagenomics humann taxonomic profile
# set $database to u50 for uniref50 analysis
# set $database to u90 for uniref90 analysis
# !! this will set the database for the whole system, not just this run !!
# !! don't try to run two sets of analysis on different databases simultaneously !!
# use search mode uniref90 for well-characterized microbiomes (ie, human gut)
# if the total reads unmapped is low or the microbiome is uncharacterized, try uniref50
# $refDir contains the paired metagenomics taxonomic profiles for metatranscriptomics analysis; can provide "" for DNA runs
# $inputDir contains the un-concatenated paired end fastq files, preferably quality trimmed with adapters removed
# $use lets the wrapper script know whether to include all the samples in $inputDir (DIRECTORY) or not (LIST)
# $list provides the list of variables to use if the LIST option is specified for $use

#### Define Run Variables
# defaults
concatenate="TRUE"
type="DNA"
refDir=""
use="DIRECTORY"
resume="FALSE"
normType="Adjusted CPMs"

VALID_ARGS=$(getopt -o i:cpo:r:n:ml:h \
                    --long inputDir:,concatenated,paired,outDir:,dnaRef:,normType:,resume,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Fastq files will be taken from the directory '$2'"
        inputDir="$2"
        shift 2
        ;;
    -c | --concatenated)
        echo "Fastq files are already concatenated"
        concatenate="FALSE"
        shift
        ;;
    -p | --paired)
        echo "Data is from a paired metatranscriptomics sample set"
        type="RNA"
        shift
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -r | --dnaRef)
        echo "Paired DNA metaphlan output is in the directory '$2'"
        refDir="$2"
        shift 2
        ;;
    -n | --normType)
        echo "Output will be normalized as '$2' values"
        normType="$2"
        shift 2
        ;;
    -m | --resume)
        echo "Will begin Humann/Metaphlan from a previous unfinished run"
        resume="TRUE"
        shift
        ;;
    -l | --list)
        echo "Samples from the list '$2' will be classified"
        list="$2"
        use="LIST"
        shift 2
        ;;
    -h | --help)
        help="TRUE"
        shift
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
  This script handles HUMAnN 4.0 analysis, including Metaphlan

  usage: sbatch C1.humann-metaphlan-wrapper.sh
            --inputDir /path/to/fastq-input
            --outDir /path/to/humann-output
            (--concatenated FALSE) (--paired FALSE)
            (--dnaRef "") (--normType "Adjusted CPMs")
            (--resume) (--list "id1 id2 id3 ...") (--help)

  options:
    [ -i  | --inputDir     | directory containing fastq.gz files                                                                              ]
    [ -c  | --concatenated | this flag should be given if the input DNA is already concatenated                                               ]
    [ -p  | --paired       | this flag should be given if the input DNA represents the RNA portion of a paired metatranscriptomics sampleset  ]
    [ -o  | --outDir       | directory where the output files will be written (including the HUMAnN temp folder)                              ]
    [ -r  | --dnaRef       | pathname to the Metaphlan profile corresponding to the DNA portion of a paired metatranscriptomics sampleset     ]
    [ -n  | --normType     | specify the normalization method for output metrics: Adjusted CPMs (default), Adjusted RPKs, RPKs, Counts        ]
    [ -m  | --resume       | this flag should be given if the run is continuing a previous run                                                ]
    [ -l  | --list         | list specific sample names (ID prior to first underscore) to run from the input directory                        ]
    [ -h  | --help         | prints an informational message and exits script                                                                 ]
EOF
  exit;
fi

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/humann4-env/

if [ "$use" == "DIRECTORY" ];
then

  ### run C1.humann-metaphlan for all samples in a directory ###

  for sid in $(find "$inputDir" -type f -name "*.fastq")
  do
    gzip "$sid"
  done

  for i in $(find "$inputDir" -type f -name "*.fastq.gz" | while read F; do basename $F; done | rev | cut -d "_" -f 5-10 | rev | sort | uniq)
  do
    echo "$i"
    sbatch --job-name "$i".C1.humann-metaphlan \
         /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_C/C1.humann-metaphlan.sh \
         "$i" "$refDir" "$outDir" "$inputDir" "$concatenate" "$type" "$mode" "$resume" "$normType"
  done;

elif [ "$use" == "LIST" ];
then

  ### to run C1.humann-metaphlan for a subset of samples in directory ###

  for i in $list
  do
    echo "$i"
    sbatch --job-name "$i".C1.humann-metaphlan \
         /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/MGX-MTX_Modules_2025-update/Module_C/C1.humann-metaphlan.sh \
         "$i" "$refDir" "$outDir" "$inputDir" "$concatenate" "$type" "$mode" "$resume" "$normType"
  done;

fi
