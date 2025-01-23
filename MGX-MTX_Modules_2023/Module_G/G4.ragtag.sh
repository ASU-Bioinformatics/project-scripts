#!/bin/bash

##### ragtag scaffolding of high-quality MAGs #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.G3.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.G3.err               # STDERR (%j = JobId)
#SBATCH -t 2-0:00                    # estimated time needed - this is very slow for these large genomes
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

module load mamba/latest

source activate /data/biocore/programs/conda-envs/ragtag-env

# create a data input file where the bin name is in the first column
# and the directory for the reference genome(s) is in the second column

binDir="/data/gencore/analysis_projects/6029563_Ball_Metagenomics/transect/3_assembly/3_MAG-bins/metabat-bins"
refDir="/data/gencore/databases/reference_genomes/rhodoferax_reference_genomes"

cd $binDir

echo $PATH

input="ragtag-test.csv"
while IFS=, read -r bin genome
do
  echo "$bin"
  echo "$genome"
  ragtag.py scaffold -o "$binDir"/"$bin"."$genome".ragtag.out \
                     "$refDir"/"$genome"/"$genome"_*.fna \
                     "$binDir"/"$bin"
done < "$input"

for i in $(find "$binDir" -type d -name "$bin.*" -maxdepth 1 | while read F; do basename $F; done | rev | cut -d '.' -f 1-3 | rev | sort | uniq)
do
  ragtag.py merge "$binDir"/"$i" "$binDir"/"$i".*.ragtag.out/*.agp other.map.agp
done;

#BUSCO quality check

#conda activate /data/biocore/programs/conda-envs/busco-env
#BUSCO_DOWNLOADS="/data/biocore/programs/busco_downloads"

#cd /data/biocore/brown-temps/busco

#busco \
#  -i /data/gencore/analysis_projects/Brown_Insect_Genome_Assembly/Liberibacter_assemblies/central/c.n-10x.libr-scaffold.fa \
#  -o c.n-10x.libr-scaffold.fa.busco.bac \
#  -m genome \
#  -l bacteria_odb10 \
#  --download_path $BUSCO_DOWNLOADS \
#  --offline

#  busco \
#    -i /data/gencore/analysis_projects/Brown_Insect_Genome_Assembly/Liberibacter_assemblies/central/c.n-10x.libr-scaffold.fa \
#    -o c.n-10x.libr-scaffold.fa.busco.proteobac \
#    -m genome \
#    -l proteobacteria_odb10 \
#    --download_path $BUSCO_DOWNLOADS \
#    --offline
