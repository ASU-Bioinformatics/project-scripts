#!/bin/bash

##### metagenome assembly using megahit #####

#SBATCH -p highmem
###SBATCH -p public
#SBATCH -q public
###SBATCH -q debug
#SBATCH --output slurm.%j.MGX.D1.out
#SBATCH --error slurm.%j.MGX.D1.err
#SBATCH -c 128
#SBATCH -t 2-0:00
###SBATCH -t 0-0:15
###SBATCH --mem=1024G

module purge
module load mamba/latest

##### Define Variables #####

resume="FALSE"
mhEnv="/data/biocore/programs/conda-envs/megahit-args"
qstEnv="/data/biocore/programs/conda-envs/quast-env"

VALID_ARGS=$(getopt -o f:r:o:p:mh \
                    --long forwardReads:,reverseReads:,outputDirectory:,outPrefix:,resume,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -f | --forwardReads)
        echo "The forward reads are: '$2'"
        forwards="$2"
        shift 2
        ;;
    -r | --reverseReads)
        echo "The reverse reads are: '$2'"
        reverses="$2"
        shift 2
        ;;
    -o | --outputDirectory)
        echo "The assembly will be built in '$2'"
        outDir="$2"
        shift 2
        ;;
    -p | --outPrefix)
        echo "The output file prefix is '$2'"
        outPrefix="$2"
        shift 2
        ;;
    -m | --resume)
        echo "MEGAHIT will continue building the assembly in '$outDir'"
        resume="TRUE"
        shift
        ;;
    -e | --megahitEnvironment)
        echo "The conda environment to use for MEGAHIT is '$2'"
        mhEnv="$2"
        shift 2
        ;;
    -q | --quastEnvironment)
        echo "The conda environment to use for QUAST is '$2'"
        qstEnv="$2"
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
  This script uses MEGAHIT to build a metagenomic assembly using paired-end short reads.
  If multiple samples should be co-assembled (if, for example, you have multiple samples from the same source),
  provide the list of forward read files to the `-forwardReads` parameter and the list of reverse read
  files to the `-reverseReads` parameter, using the same sample order for each.

  usage: sbatch D1.megahit-assemble.sh
            -f "A.r1.fastq.gz B.r1.fastq.gz" OR -f "forwards.txt"
            -r "A.r2.fastq.gz B.r2.fastq.gz" OR -f "reverses.txt"
            (-o /path/to/output/directory) (-p output.prefix)
            (-e /path/to/megahit/environment)
            (-q /path/to/quast/environment) (-m) (-h)

  options:
    [ -f  |   --forwardReads        | list of, or text file containing, the pathnames for all forward read files              ]
    [ -r  |   --reverseReads        | list of, or text file containing, the pathnames for all reverse read files              ]
    [ -o  |   --outputDirectory     | subfolder where fastq files will be moved following QC (default 'fastq')                ]
    [ -p  |   --outPrefix           | prefix for the final contigs file                                                       ]
    [ -m  |   --resume              | if provided, the script will resume an unfinished assembly                              ]
    [ -e  |   --megahitEnvironment  | megahit environment location; default "/data/biocore/programs/conda-envs/megahit-args"  ]
    [ -q  |   --quastEnvironment    | quast environment location; default "/data/biocore/programs/conda-envs/quast-env"       ]
    [ -h  |   --help                | prints an informational message and exits script                                        ]
EOF
  exit;
fi

#### run megahit #####

source activate "$mhEnv"

if [[ "$resume" == "FALSE" ]]; then
  megahit -1 $forwards \
          -2 $reverses \
          -o "$outDir" \
          --presets meta-large -t 128 \
          --out-prefix $outPrefix
else
  megahit -1 $forwards \
          -2 $reverses \
          -o "$outDir" --continue \
          --presets meta-large -t 128 \
          --out-prefix $outPrefix
fi

source deactivate
source activate "$qstEnv"

cd "$outDir"

quast -o "$outPrefix"-dna-quast "$outDir"/"$outPrefix".contigs.fa
