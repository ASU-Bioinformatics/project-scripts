#!/bin/bash

#SBATCH -o slurm.%j.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err               # STDERR (%j = JobId)

### for sol ###
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-01:00:00
#SBATCH --mem=16G

module load mamba/latest

source activate /data/biocore/programs/conda-envs/trapp-env

VALID_ARGS=$(getopt -o i:o: \
                    --long inputDir:,outDir: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -i | --inputDir)
        echo "Will look in the directory '$2' for raw Bracken files"
        inputDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
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

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script uses the alpha-diversity python script from KrakenTools 1.2 to analyze all Bracken reports contained in the input directory.

  usage: sbatch B3.bracken-diversity.sh
            --inputDir /path/to/bracken-reports
            --outDir /path/to/bracken/diversity-output
            (-h)

  options:
    [ -i  |   --inputDir  |   directory containing the kraken report files (as output by B1.kraken.sh, file name ending in .k2report)  ]
    [ -o  |   --outDir    |   directory for bracken output folders (raw output and reports)                                            ]
    [ -h  |   --help      |   prints an informational message and exits script                                                         ]
EOF
  exit;
fi

mkdir -p "$outDir"

cd "$inputDir"

echo "Shannon
      Berger-Parker
      Simpson
      InverseSimpson
      Fisher" > "$outDir"/summary.txt

for i in $(find ./ -type f -name "*.bracken" | while read F; do basename $F; done)
do
  echo "$i" >> "$i".txt
  /data/biocore/programs/KrakenTools-1.2/DiversityTools/alpha_diversity.py \
    -f "$i" -a Sh | grep '[[:digit:]]' | grep -o '[0-9.]*' >> "$outDir"/"$i".txt
  /data/biocore/programs/KrakenTools-1.2/DiversityTools/alpha_diversity.py \
    -f "$i" -a BP | grep '[[:digit:]]' | grep -o '[0-9.]*' >> "$outDir"/"$i".txt
  /data/biocore/programs/KrakenTools-1.2/DiversityTools/alpha_diversity.py \
    -f "$i" -a Si | grep '[[:digit:]]' | grep -o '[0-9.]*' >> "$outDir"/"$i".txt
  /data/biocore/programs/KrakenTools-1.2/DiversityTools/alpha_diversity.py \
    -f "$i" -a ISi | grep '[[:digit:]]' | grep -o '[0-9.]*' >> "$outDir"/"$i".txt
  /data/biocore/programs/KrakenTools-1.2/DiversityTools/alpha_diversity.py \
    -f "$i" -a F | grep '[[:digit:]]' | grep -o '[0-9.]*' >> "$outDir"/"$i".txt
done

cd "$outDir"

txt=$(find ./ -type f -name "*S.bracken.txt" | while read F; do basename $F; done | sort | uniq)
paste summary.txt $txt > "$outDir"/species_alpha-diversity.txt

txt=$(find ./ -type f -name "*G.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/genus_alpha-diversity.txt

txt=$(find ./ -type f -name "*F.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/family_alpha-diversity.txt

txt=$(find ./ -type f -name "*O.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/order_alpha-diversity.txt

txt=$(find ./ -type f -name "*C.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/class_alpha-diversity.txt

txt=$(find ./ -type f -name "*P.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/phylum_alpha-diversity.txt

txt=$(find ./ -type f -name "*K.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/kingdom_alpha-diversity.txt

txt=$(find ./ -type f -name "*D.bracken.txt" | while read F; do basename $F; done| sort | uniq)
paste summary.txt $txt > "$outDir"/domain_alpha-diversity.txt

cd "$inputDir"

beta=$(find ./ -type f -name "*S.bracken" | while read F; do basename $F; done | sort | uniq)

# creates bray-curtis dissimilarity matrix - 0=identical composition, 1=no species shared
/data/biocore/programs/KrakenTools-1.2/DiversityTools/beta_diversity.py -i $beta --type bracken --level all > beta-diversity.txt
mv beta-diversity.txt "$outDir"
