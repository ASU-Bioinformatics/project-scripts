#!/bin/bash

##### merge microbe-annotator runs from split predicted protein files #####

#SBATCH -p htc
#SBATCH -q public
#SBATCH -o slurm.%j.F4.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.F4.err               # STDERR (%j = JobId)
#SBATCH -t 0-04:00
#SBATCH -c 4
#SBATCH --mem=32G

echo "merge microbe annotator final annot and ko files"

outDir=$1

for i in $(find "$outDir" -type f -name "out*.annot" | while read F; do basename "$F"; done)
do
  echo "$i"
  tail -n +2 "$outDir"/microbe-annotator-out-*/annotation_results/"$i" > "$outDir"/noheader-"$i"
done

for i in $(find "$outDir" -type f -name "noheader*annot");
do
  echo "$i"
  cat "$i" >> "$outDir"/all.out.CAT.predicted_proteins.faa.annot
done;

header=$(head -n 1 "$outDir"/microbe-annotator-out-00/annotation_results/out.CAT.predicted_proteins.00.faa.annot)
echo "$header" > "$outDir"/header-annot.txt
cat "$outDir"/header-annot.txt "$outDir"/all.out.CAT.predicted_proteins.faa.annot > "$outDir"/final.out.CAT.predicted_proteins.faa.annot

for i in $(find "$outDir" -type f -name "*.ko" | while read F; do basename "$F"; done)
do
  echo "$i"
  tail -n +2 "$outDir"/microbe-annotator-out-*/annotation_results/"$i" > "$outDir"/noheader-"$i"
done

for i in $(find "$outDir" -type f -name "noheader*ko");
do
  echo "$i"
  cat "$i" >> "$outDir"/all.out.CAT.predicted_proteins.faa.ko
done;

header=$(head -n 1 "$outDir"/microbe-annotator-out-00/annotation_results/out.CAT.predicted_proteins.00.faa.ko)
echo "$header" > "$outDir"/header-ko.txt
cat "$outDir"/header-ko.txt "$outDir"/all.out.CAT.predicted_proteins.faa.ko > "$outDir"/final.out.CAT.predicted_proteins.faa.ko

mkdir -p "$outDir"/unmerged-annotations
mv microbe-annotator-out-* "$outDir"/unmerged-annotations/
mv "$outDir"/noheader-out.CAT.predicted_proteins.0* "$outDir"/unmerged-annotations/
mv out.CAT.predicted_proteins.*.faa "$outDir"/unmerged-annotations/
mv header* "$outDir"/unmerged-annotations/
