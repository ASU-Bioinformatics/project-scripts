#!/bin/bash

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%A.out               # STDOUT (%A = JobId, %a = TaskID)
#SBATCH -e slurm.%A.err               # STDERR (%A = JobId, %a = TaskID)
#SBATCH -t 0-4:00
#SBATCH -c 1
#SBATCH --mem=16G

module purge
module load mamba/latest

source activate /data/biocore/programs/mamba-envs/biocore-rna
PATH=/data/biocore/programs/mamba-envs/biocore-rna/bin/:$PATH

help="FALSE"
useList="FALSE"

VALID_ARGS=$(getopt -o i:d:l:h \
                    --long infoFile:,inputDirectory:,list:,help \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --infoFile)
        echo "the gene information for annotation is in '$2'"
        infoFile="$2"
        shift 2
        ;;
    -d | --inputDirectory)
        echo "Files for annotation are located in '$2'"
        inputDirectory="$2"
        shift 2
        ;;
    -l | --list)
        echo "Only samples in the list '$2' within the input directory will be annotated"
        list=$2
        useList="TRUE"
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
  This script combines gene information with ID'd-only gene lists to annotate the gene lists.
  This improves readability and usefulness for downstream research.

  usage: bash C2.annotate-DEGs \
            --infoFile /path/to/geneInfo.tab \
            --inputDirectory /path/to/deglists \
            (--list "file1 file2 file3")

  options:
    [ -i  |   --infoFile        |   file containing gene IDs matching the DEG list IDs along with any desired information  ]
    [ -d  |   --inputDirectory  |   directory containing the files                                                         ]
    [ -l  |   --list            |   quoted, space-separated list of files within the input directory to annotate           ]
    [ -h  |   --help            |   prints informational message and exits                                                 ]

EOF
exit;
fi

# determine list of files to annotate
# need to find out how to exclude image files when extensions are unknown and/or variable
if [ "$useList" == "FALSE" ]; then
  list=$(find "$inputDirectory" -maxdepth 1 -type f -name "*" | while read F; do basename $F; done)
fi

# find delimiter for infoFile
# using tail helps avoid any annotation or header row that might cause problems
# for example, Star index gene information files have an extra number at the top of the file.
inum=1

for delim in '\t' ',' ';' ' ' '|'; do
  echo "$delim"
  inum=$(tail -n 1 $infoFile | tr "$delim" '\n' | wc -l )
  echo "$inum"
  if [ $inum -gt 1 ]; then
    echo "delimiter of $infoFile is $delim"
    break;
  else
    continue;
  fi
done

cd "$inputDirectory"

infoName=$(basename $infoFile)
tr "$delim" '\t' < $infoFile > delim."$infoName"
wc -l $infoFile
wc -l delim."$infoName"
sortInfo=sorted.delim."$infoName"
LC_ALL=C
LANG=C
sort -t $'\t' delim."$infoName" > "$sortInfo"

num=1
found="FALSE"
for file in $list; do
  head $file
  for delim in '\t' ',' ';' ' ' '|'; do
    echo "$delim"
    num=$(tail -n 1 $file | tr "$delim" '\n' | wc -l )
    echo "$num"
    if [ $num -gt 1 ]; then
      echo "delimiter of $file is $delim"
      found="TRUE"
      break;
    fi
  done
  if [[ "$found" == "FALSE" ]]; then
    exit;
  fi
  tr -d '"' < $file | sed -z 's/\n'$delim'/'$delim'/g' | tr "$delim" '\t' > delim."$file"
  #tr "$delim" '\t' < $file > delim."$file"
  head delim."$file"
  header=""
  for ((i=0; i<"$inum"; i++)); do
    header+="null"$'\t';
  done
  echo $header
  header+=$(head -n 1 delim."$file")
  echo $header
  annFile=ann."$file".txt
  LC_ALL=C
  LANG=C
  sort -t $'\t' delim.$file > sorted.delim.$file
  wc -l sorted.delim.$file
  head sorted.delim.$file
  wc -l $sortInfo
  #head $sortInfo
  join -t $'\t' $sortInfo sorted.delim.$file > $annFile
  head $annFile
  sed -i "1i $header" $annFile
done
