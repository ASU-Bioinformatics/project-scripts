#!/bin/bash

##### humann3 test #####
#### %j is job id, %x is job name (specified in loop to call script)

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH -t 3-12:00
#SBATCH -c 6
#SBATCH --mem=84G

sid="$1"

echo "humann run for $sid, set for 3.5 days, 6 cores, 84G memory"

module load mamba/latest
source activate /data/biocore/programs/conda-envs/humann-env/

refDir="$2"
outDir="$3"
fqDir="$4"
concatenate="$5"
type="$6"
database="$7"

if [[ "$concatenate" == "TRUE" ]]; then
  echo "concatenating fastq files"
  mkdir -p "$outDir"
  mkdir -p "$fqDir"/concatenated
  cat "$fqDir"/"$sid"_SQP_L00*_R1_001.fastq.gz "$fqDir"/"$sid"_SQP_L00*_R2_001.fastq.gz > "$fqDir"/concatenated/"$sid"_SQP_L001_RC_001.fastq.gz
  fqDir="$fqDir"/concatenated
fi

fastq="$fqDir"/"$sid"_SQP_L001_RC_001.fastq.gz

if [[ "$database" == "u50" ]]; then
  humann_config --update database_folders protein /data/gencore/databases/humann/uniref/uniref
  mode="uniref50"
elif [[ "$database" == "u90" ]]; then
  humann_config --update database_folders protein /data/gencore/databases/humann/uniref90/
  mode="uniref90"
fi

if [[ "$type" == "DNA" ]]; then
  echo "running DNA or unpaired RNA mode"
  humann -i "$fastq" --search-mode "$mode" \
         -o "$outDir" --output-format tsv \
         --threads 8 --verbose
elif [[ "$type" == "RNA" ]]; then
  echo "running paired RNA mode"
  humann -i "$fastq" --search-mode "$mode" \
         -o "$outDir" --output-format tsv \
         --threads 8 --verbose \
         --taxonomic-profile "$refDir"/"$sid"_SQP_L001_RC_001_humann_temp/"$sid"_SQP_L001_RC_001_metaphlan_bugs_list.tsv
fi


#take metaphlan output in tmp folders and run MaAsLin analysis and make a heatmap of top variable species
#can do the same with humann3 manipulated output from C2 script for gene families and pathway abundances
