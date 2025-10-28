#!/bin/bash
#SBATCH -c 1
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-4:00
#SBATCH --mem=100G
#SBATCH -o slurm.%A.markdup.out
#SBATCH -e slurm.%A.markdup.err

# markdup takes more memory than I easily can have in an interactive session
# and I'm not sure how long it will take
# I'm going to try the first vaccinia with 4hrs and 100G and see how that goes.

umask 0007

module load samtools-1.21-gcc-12.1.0

mDir="/data/gencore/analysis_projects/8406252_Bonner/fastq/mouse"
vDir="/data/gencore/analysis_projects/8406252_Bonner/fastq/vaccinia"

#samtools markdup -r -f $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
#                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
#                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.bam

#samtools markdup -r -f $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
#                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
#                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.bam

#samtools markdup -r -f $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.stats.txt \
#                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
#                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.bam

#samtools markdup -r -f $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.stats.txt \
#                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
#                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam

#samtools markdup -r -f $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.stats.txt \
#                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
#                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.bam

samtools markdup -r -f $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.stats.txt \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam
