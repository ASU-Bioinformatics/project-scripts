#!/bin/bash

##### trying to get abundance tables #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.%x.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.%x.err               # STDERR (%j = JobId)
#SBATCH -t 0-8:00
#SBATCH -c 1
#SBATCH --mem=64G

module load mamba/latest

source deactivate
source activate /data/biocore/programs/mamba-envs/htseq-env

VALID_ARGS=$(getopt -o a:o:g:l: \
                    --long alignmentDir:,outDir:,gffFile:,list: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -a | --alignmentDir)
        echo "Transcript alignments will be read from '$2'"
        alnDir="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "Results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -g | --gffFile)
        echo "Coassembly annotations will be read from '$2'"
        gffFile="$2"
        shift 2
        ;;
    -l | --list)
        echo "Samples will be taken from the following list of IDs present in the transcript directory: '$2'"
        use=LIST
        list="$2"
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

sid="$1"

echo "$sid htseq quantification using CAT predicted proteins, 8 hours, 1024G memory, 16 threads"

alnDir="$2"

gffFile="$3"

outDir="$4"

htseq-count -f bam -r pos -s no -t CDS -i ID --add-chromosome-info \
            -d "\t" -n 16 --max-reads-in-buffer=30000000000000 \
            "$alnDir"/"$sid".sorted.bam \
            "$gffFile" > "$outDir"/"$sid".htseqcounts.txt

source deactivate
module load samtools-1.13-gcc-11.2.0

samtools flagstat "$alnDir"/"$sid".headed.bam > "$sid".flagstats.txt
