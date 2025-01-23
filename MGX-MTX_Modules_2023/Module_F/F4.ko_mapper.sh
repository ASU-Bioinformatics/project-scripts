#!/bin/bash

##### download and test interproscan #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.F4.out
#SBATCH -o slurm.%j.F4.err
#SBATCH -t 0-4:00
#SBATCH --mem=64G

#### wrap-up low memory scripts for CAT annotations ####

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/catbat-env

#### official default values ####
taxDir="./CAT_taxonomy"
scriptDir="./"
contigs="contigs.fna"
outDir="./CAT_output"
help="FALSE"

#### default variables for core use (just so I don't forget them); provide these through the command line ####
#taxDir="/data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy"
#scriptDir="/data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline"

VALID_ARGS=$(getopt -o t:s:c:o:h \
                    --long taxDir:,scriptDir:,contigs:,outDir:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -t | --taxDir)
        echo "CAT taxonomy files are located in '$2'."
        taxDir="$2"
        shift 2
        ;;
    -s | --scriptDir)
        echo "MicrobeAnnotator's ko_mapper.py script is located in '$2'."
        scriptDir="$2"
        shift 2
        ;;
    -c | --contigs)
        echo "The contig assembly file is '$2'"
        contigs="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Output files will be written to '$2'"
        outDir="$2"
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

if [[ $help == "TRUE" ]];

then
  echo "runs several auxiliary scripts to summarize CAT and MicrobeAnnotator outputs"
  echo "takes as input the out.CAT.contig2classification.txt file from CAT, the final assembly, and the CAT taxonomy directory"
  echo "usage showing default values for optional arguments:"
  echo "sbatch F4.ko_mapper.sh -t ./CAT_taxonomy -s ./ -c contigs.fa -o ./CAT_output"
  echo "all samples from the alignment folder will be used unless the -l option provides a quoted, space-separated list of sample ids"

else
#### CAT classify with official names ####

  mkdir -p "$outDir"

  CAT add_names -i "$outDir"/out.CAT.contig2classification.txt \
                -o "$outDir"/out.CAT.contig2classification.official.txt \
                -t "$taxDir" --only_official

#### CAT summarise ####

  CAT summarise -c "$contigs" \
                -i "$outDir"/out.CAT.contig2classification.official.txt \
                -o "$outDir"/out.CAT.summary.txt

#### map list of Kegg terms from CAT with ko_mapper from MicrobeAnnotator ####

  source deactivate

  source activate /data/biocore/programs/mamba-envs/microbe-annotator-env

  cd "$outDir"

  python "$scriptDir"/ko_mapper.py --input_files "$outDir"/final.out.CAT.predicted_proteins.faa.ko \
                                   --prefix CAT.contigs --cluster rows
fi
