#!/bin/bash

##### humann3 test #####
#### %j is job id, %x is job name (specified in loop to call script)

#SBATCH -p public
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out
#SBATCH -e slurm.%j.%x.err
#SBATCH -t 7-00:00
#SBATCH -c 4
#SBATCH --mem=480G

sid="$1"

echo "humann run for $sid, set for 3.5 days, 6 cores, 84G memory"

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/humann4-env/

refDir="$2"
outDir="$3"
fqDir="$4"
concatenate="$5"
type="$6"
mode="$7"
resume="$8"

if [[ "$concatenate" == "TRUE" ]]; then
  echo "concatenating fastq files"
  mkdir -p "$outDir"
  mkdir -p "$fqDir"/concatenated
  cat "$fqDir"/"$sid"_SQP_L00*_R1_001.fastq.gz "$fqDir"/"$sid"_SQP_L00*_R2_001.fastq.gz > "$fqDir"/concatenated/"$sid"_SQP_L001_RC_001.fastq.gz
  fqDir="$fqDir"/concatenated
fi

fastq="$fqDir"/"$sid"_SQP_L001_RC_001.fastq.gz


if [[ "$type" == "DNA" ]]; then
  echo "running DNA or unpaired RNA mode"

  if [[ "$resume" == "TRUE" ]]; then
    echo "resuming from previous run"
    humann -i "$fastq" --metaphlan-options="--offline -t rel_ab_w_read_stats" \
           -o "$outDir" --output-format tsv \
           --threads 6 --verbose --resume
  else
    humann -i "$fastq" --metaphlan-options="--offline -t rel_ab_w_read_stats" \
          -o "$outDir" --output-format tsv \
          --threads 6 --verbose
  fi

elif [[ "$type" == "RNA" ]]; then
  echo "running paired RNA mode"

  if [[ "$resume" == "TRUE" ]]; then
    echo "resuming from previous run"
    humann -i "$fastq" --metaphlan-options="--offline -t rel_ab_w_read_stats" \
           -o "$outDir" --output-format tsv \
           --threads 6 --verbose \
           --taxonomic-profile "$refDir"/"$sid"_SQP_L001_RC_001_humann_temp/"$sid"_SQP_L001_RC_001_metaphlan_bugs_list.tsv
  else
    humann -i "$fastq" --metaphlan-options="--offline -t rel_ab_w_read_stats" \
          -o "$outDir" --output-format tsv \
          --threads 6 --verbose \
          --taxonomic-profile "$refDir"/"$sid"_SQP_L001_RC_001_humann_temp/"$sid"_SQP_L001_RC_001_metaphlan_bugs_list.tsv
  fi
fi


#take metaphlan output in tmp folders and run MaAsLin analysis and make a heatmap of top variable species
#can do the same with humann3 manipulated output from C2 script for gene families and pathway abundances
