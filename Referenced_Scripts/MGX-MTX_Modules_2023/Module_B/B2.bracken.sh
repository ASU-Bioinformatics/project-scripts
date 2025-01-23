#!/bin/bash

#SBATCH -o slurm.%j.B2.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.B2.err               # STDERR (%j = JobId)

### for sol ###
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-04:00
#SBATCH --mem=4G

levels="S G F O C P K D"

VALID_ARGS=$(getopt -o i:o:d: \
                    --long inputDir:,outDir:,dbDir: \
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
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -d | --dbDir)
        echo "The KRAKEN database in the directory '$2' will be used for classification"
        use="$2"
        shift 2
        ;;
    -t | --levels)
        echo "Bracken abundance will be calculated for levels '$2'"
        levels="$2"
        shift 2
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

dbDir="/data/gencore/databases/kraken/k2_pluspf"

brkDir="$outDir/bracken-raws"
brpDir="$outDir/bracken-reports"

mkdir -p "$brkDir"
mkdir -p "$brpDir"

for level in $levels;
do
  echo "$level"
  for i in $(find "$inputDir" -maxdepth 1 -type f -name "*.k2report" | while read F; do basename $F; done | sort | uniq);
  do
    echo "$i"
    echo ${i%.k2report}_"$level".breport
    /data/biocore/programs/bracken/Bracken-2.8/bracken -d "$dbDir" -i "$inputDir"/"$i" -r 150 -l "$level" -t 10 \
    -o "$brkDir"/${i%.k2report}_"$level".bracken -w "$brpDir"/${i%.k2report}_"$level".breport
  done;
done;
