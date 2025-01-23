#!/bin/bash

##### ARG-focused metagenomics analysis using megahit #####

#SBATCH -p highmem
#SBATCH -q public
#SBATCH -o slurm.%j.D2.out              # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.D2.err              # STDERR (%j = JobId)
#SBATCH -c 128
#SBATCH -t 2-0:00
#SBATCH --mem=1024G

module purge
module load mamba/latest

##### Define Variables #####

assemblyName="final.contigs.fa"

VALID_ARGS=$(getopt -o a:p:n: \
                    --long assemblyDir:,prefix:,assemblyName: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --assemblyDir)
        echo "The assembly located in '$2' will be indexed."
        assemblyDir="$2"
        shift 2
        ;;
    -p | --prefix)
        echo "The prefix for the bowtie2 index will be '$2'"
        prefix="$2"
        shift 2
        ;;
    -n | --assemblyName)
        assemblyName="$2"
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

echo "The assembly name in '$assemblyDir' is '$assemblyName'"

#### build bowtie2 index ####

assembly="$assemblyDir"/"$assemblyName"

source activate /data/biocore/programs/conda-envs/bowtie2-env

bowtie2-build --threads 64 "$assembly" "$assemblyDir"/"$prefix"

source deactivate

#### create assembly graphs with gfastats ####

source activate /data/biocore/programs/mamba-envs/gfastats-env

gfaDir="$assemblyDir"/assembly_graphs
mkdir "$gfaDir"

gfastats "$assemblyDir"/final.contigs.fa -o gfa > "$gfaDir"/final.contigs.gfastats.gfa
gfastats "$gfaDir"/final.contigs.gfastats.gfa --nstar-report > "$gfaDir"/final.contigs.nstarreport.txt
gfastats "$gfaDir"/final.contigs.gfastats.gfa --seq-report > "$gfaDir"/final.contigs.seqreport.txt

source deactivate
