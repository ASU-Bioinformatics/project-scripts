#!/bin/bash

##### identifying viral contigs with VIBRANT #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.E1.vibrant.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.E1.vibrant.err               # STDERR (%j = JobId)
#SBATCH -t 3-00:00
#SBATCH -c 32
#SBATCH --mem=120G

#### this script runs VIBRANT on a contig assembly file to predict and extract viral contigs and perform community functional analysis
# outDir is the absolute pathname to the directory where output files will be written
# contigs is the absolute pathname to the assembly fasta file
# databaseDir is the absolute pathname to the databases subfolder within the VIBRANT directory
# modelDir is the absolute pathname to the model subfolder within the VIBRANT directory

#### set up environment ####

module load mamba/latest

source activate /data/biocore/programs/mamba-envs/vibrant-env

#### read in arguments from command line ####

VALID_ARGS=$(getopt -o c:d:m:o: \
                    --long contigs:,databaseDir,modelDir,outDir: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -c | --contigs)
        echo "Contig file for viral prediction is '$2'"
        contigs="$2"
        shift 2
        ;;
    -d | --databaseDir)
        echo "VIBRANT database files are in '$2'"
        databaseDir="$2"
        shift 2
        ;;
    -m | --modelDir)
        echo "VIBRANT model files are in '$2'"
        modelDir="$2"
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


#### download databases (only do once) ####

#cd /data/gencore/databases

#download-db.sh /data/gencore/databases/vibrant2

# this script did not compress the hmm files correctly
# following the help on GitHub, I deleted the .h3i files and ran hmmpress on each of the HMM files

#hmmpress KEGG_profiles_prokaryotes.HMM
#hmmpress Pfam-A_v32.HMM
#hmmpress VOGDB94_phage.HMM

#### run vibrant on contigs file ####

prefix=$(basename $contigs | rev | cut -d '.' -f 2-10 | rev)

mkdir -p "$outDir"
cd "$outDir"

VIBRANT_run.py -t 128 -f nucl -folder "$outDir" \
               -virome -d "$databaseDir" -m "$modelDir" \
               -i "$contigs"

prodigal -f gff -o "$outDir"/"$prefix".prodigal.gff \
         -i "$outDir"/VIBRANT_final.contigs/VIBRANT_phages_"$prefix"/"$prefix".phages_combined.fna \
         -a "$outDir"/VIBRANT_final.contigs/VIBRANT_phages_"$prefix"/"$prefix".phages_combined.faa \
         -d "$outDir"/VIBRANT_final.contigs/VIBRANT_phages_"$prefix"/"$prefix".phages_combined.ffn \
         -p meta
