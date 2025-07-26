# clip-seq workflow

## download reference genomes

These genomes are the primary mouse reference genome GRCm39 from Ensembl and the species exemplar for *Vaccinia*, accessed via NCBI.

```
wget https://ftp-ncbi-nlm-nih-gov.ezproxy1.lib.asu.edu/genomes/all/GCF/000/860/085/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.gtf.gz
wget https://ftp-ncbi-nlm-nih-gov.ezproxy1.lib.asu.edu/genomes/all/GCF/000/860/085/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.gff.gz
wget https://ftp-ncbi-nlm-nih-gov.ezproxy1.lib.asu.edu/genomes/all/GCF/000/860/085/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.gbff.gz
wget https://ftp-ncbi-nlm-nih-gov.ezproxy1.lib.asu.edu/genomes/all/GCF/000/860/085/GCF_000860085.1_ViralProj15241/GCF_000860085.1_ViralProj15241_genomic.fna.gz

wget https://ftp.ensembl.org/pub/release-114/gff3/mus_musculus/Mus_musculus.GRCm39.114.gff3.gz
wget https://ftp.ensembl.org/pub/release-114/gtf/mus_musculus/Mus_musculus.GRCm39.114.gtf.gz
wget https://ftp.ensembl.org/pub/release-114/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz
```

## data pre-processing

### adaptor trimming with cutadapt

This initial script is taken from the PureCLIP documentation. Since our libraries weren't built with a CLIP-specific protocol, I am not sure whether it would be correct to use these sequences or if I should use the standard adaptor file I use for our in-house Illumina libraries.

```
cutadapt -f fastq --match-read-wildcards --times 1 -e 0.1 -O 1 \
         --quality-cutoff 6 -m 18 -a NNNNNAGATCGGAAGAGCACACGTCTGAACTCCAGTCAC \
         -g CTTCCGATCTACAAGTT  -g CTTCCGATCTTGGTCCT \
         -A AACTTGTAGATCGGA -A AGGACCAAGATCGGA -A ACTTGTAGATCGGAA -A GGACCAAGATCGGAA  \
         -A CTTGTAGATCGGAAG  -A GACCAAGATCGGAAG -A TTGTAGATCGGAAGA -A ACCAAGATCGGAAGA \
         -A TGTAGATCGGAAGAG -A CCAAGATCGGAAGAG -A GTAGATCGGAAGAGC -A CAAGATCGGAAGAGC \
         -A TAGATCGGAAGAGCG  -A AAGATCGGAAGAGCG -A AGATCGGAAGAGCGT  -A GATCGGAAGAGCGTC \
         -A ATCGGAAGAGCGTCG -A TCGGAAGAGCGTCGT -A CGGAAGAGCGTCGTG -A GGAAGAGCGTCGTGT \
         -o rep1/reads.R1.trimmed.fastq.gz -p rep1/reads.R2.trimmed.fastq.gz \
         rep1/reads.R1.fastq.gz rep1/reads.R2.fastq.gz
```

If I use our in-house adaptors, while maintaining as many of the above parameters as possible, the code will look like the following code block. I think this option will be correct, and it is what I will try first.

To make the sample names compatible with downstream code, I renamed them as follows:

`rename 'L008' 'SRN_L001' *`

```
#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.cutadapt.out
#SBATCH -e slurm.%A.cutadapt.err
#SBATCH -t 0-12:00
#SBATCH --mem=32G

umask 0007

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/cutadapt/

fDir=/data/gencore/analysis_projects/8406252_Bonner/fastq

cd $fDir

adapters="/data/gencore/databases/trimmomatic/nextera.fa"

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  cutadapt -a file:"$adapters" -A file:"$adapters" \
           --match-read-wildcards --times 1 \
           -e 0.1 -O 1 --quality-cutoff 6 -m 18 \
           -o "$i"_SCT_L001_R1_001.fastq.gz \
           -p "$i"_SCT_L001_R2_001.fastq.gz \
           "$i"_S*_L001_R1_001.fastq.gz "$i"_S*_L001_R2_001.fastq.gz
done
```

### read mapping with star

This step aligns the filtered reads to a reference genome. In the case of a host-pathogen system, it can align to a joint genome instead of two separate reference genomes if the joint genome index is generated in STAR.

Since STAR requires the reference fasta file and a reference GTF file, creating a joint genome requires concatenating the two reference fasta files into one and the two reference GTF files into one. However, it is possible to run alignment separately, which is what I plan to do since PureCLIP recommends using end-to-end mapping and eliminating multi-mappers.

#### creation of star genome indexes

I do plan to seff these jobs after they complete to get a better idea of what time, memory, and cores parameters are more suited to generating STAR indexes for both large genomes such as mouse and small viral genomes such as *Vaccinia*.

```
#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only as phx doesn't have an updated STAR module.
#SBATCH -t 0-8:00
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.mouse.out
#SBATCH -o slurm.%j.starGG.mouse.err

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/mouse/GRCm39"
fasta="$refDir"/"Mus_musculus.GRCm39.dna.primary_assembly.fa"
gtf="$refDir"/"Mus_musculus.GRCm39.114.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon exon \
  --limitGenomeGenerateRAM 3000000000000
```

For *Mus musculus*, the CPU utilized was only 02:11:26, using only 8.32% of the 1-02:20:24 core walltime (in other words, only 4 cores vs the requested 36 is needed, and even that might be overkill). The wall-clock time was only 00:43:54, so only 1hr needed to be requested. Finally, the memory utilized was 76.46G, less than the 100G requested, though the 100G may still be recommended in case memory usage for a similar-sized genome runs higher than this.

```
#!/bin/bash
#SBATCH -c 36
#SBATCH -N 1
#SBATCH -p general
#SBATCH -q public #sol only, phx doesn't have an updated STAR module.
#SBATCH -t 0-8
#SBATCH --mem=100G
#SBATCH -o slurm.%j.starGG.vaccinia.out
#SBATCH -o slurm.%j.starGG.vaccinia.err


module load mamba/latest
source activate /data/biocore/programs/mamba-envs/biocore-rna

refDir="/data/gencore/databases/reference_genomes/vaccinia"
fasta="$refDir"/"GCF_000860085.1_ViralProj15241_genomic.fna"
gtf="$refDir"/"GCF_000860085.1_ViralProj15241_genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon CDS \
  --limitGenomeGenerateRAM 3000000000000
```

For *Vaccinia*, the CPU utilized was only 00:00:06, using only 0.69% of the 00:14:24 core walltime (in other words, only 1 core vs the requested 36 is needed). The wall-clock time was only 00:00:24, so only 1hr needed to be requested. Finally, the memory utilized was only 2.94G, far less than the 100G requested.

#### alignment to reference genomes

The code here has some parameters specifically optimized for CLIP sequencing - hopefully a lot of reads pass through these more stringent filters!

This first code block contains the basic template for the sbatch code. I ran each sample against both reference genomes.

```
#!/bin/bash

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%A.align.out
#SBATCH -e slurm.%A.align.err
#SBATCH -t 0-12:00
#SBATCH -c 6
#SBATCH --mem=32G

id=$1
refDir=$2
speciesID=$3

umask 0007
module purge

ulimit -v 30000000000

fastqDir="/data/gencore/analysis_projects/8406252_Bonner/fastq"
intDir=/scratch/kawoodbu/Bonner-"$id"-"$speciesID"
read1=$fastqDir/"$id"_SCT_L001_R1_001.fastq.gz
read2=$fastqDir/"$id"_SCT_L001_R2_001.fastq.gz
prefix="$id"_STAR

mkdir -p $intDir

cd $fastqDir

cd "$intDir"

STAR --outSAMtype BAM SortedByCoordinate \
     --runThreadN 10 --genomeDir $refDir \
     --readFilesIn $read1 $read2 \
     --readFilesCommand gunzip -c --outFilterType BySJout --outFilterMultimapNmax 1 \
     --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 \
     --outFilterMismatchNoverLmax 0.04 --scoreDelOpen -1 \
     --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 \
     --outFileNamePrefix $prefix --alignEndsType EndToEnd

mkdir -p $fastqDir/"$speciesID"

mv * $fastqDir/"$speciesID"
```

```
sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh 37N-FLAG-IP /data/gencore/databases/reference_genomes/mouse/GRCm39 mouse
sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh 37N-FLAG-IP /data/gencore/databases/reference_genomes/vaccinia vaccinia

sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh Mock-FLAG-IP /data/gencore/databases/reference_genomes/mouse/GRCm39 mouse
sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh Mock-FLAG-IP /data/gencore/databases/reference_genomes/vaccinia vaccinia

sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh 37N-Total-DI /data/gencore/databases/reference_genomes/mouse/GRCm39 mouse
sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh 37N-Total-DI /data/gencore/databases/reference_genomes/vaccinia vaccinia

sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh Mock-Total-DI /data/gencore/databases/reference_genomes/mouse/GRCm39 mouse
sbatch /data/gencore/shared_scripts/RNAseq/CLIP-RNASeq/A2.align.sh Mock-Total-DI /data/gencore/databases/reference_genomes/vaccinia vaccinia
```

### filter alignment files

#### keep main chromosome reads only

First, filter the reads to get a bam file containing only the reads mapped to the main chromosomes. For mouse, these are chromosomes 1-19, X, and Y (we're excluding the mitochrondrial genes also).

```
module load samtools-1.21-gcc-12.1.0
mDir="/data/gencore/analysis_projects/8406252_Bonner/fastq/mouse"
vDir="/data/gencore/analysis_projects/8406252_Bonner/fastq/vaccinia"

# mouse
samtools index $mDir/Mock-FLAG-IP-mouse_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $mDir/Mock-FLAG-IP-mouse_STARAligned.sortedByCoord.out.bam \
              -o $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.bam \
              1:1 2:1 3:1 4:1 5:1 6:1 7:1 8:1 9:1 10:1 11:1 12:1 13:1 14:1 15:1 16:1 17:1 18:1 19:1 X:1 Y:1
samtools index $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.bam

samtools index $mDir/37N-Total-DI-mouse_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $mDir/37N-Total-DI-mouse_STARAligned.sortedByCoord.out.bam \
              -o $mDir/37N-Total-DI-mouse_STARAligned.filtered.bam \
              1:1 2:1 3:1 4:1 5:1 6:1 7:1 8:1 9:1 10:1 11:1 12:1 13:1 14:1 15:1 16:1 17:1 18:1 19:1 X:1 Y:1
samtools index $mDir/37N-Total-DI-mouse_STARAligned.filtered.bam

samtools index $mDir/37N-FLAG-IP-mouse_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $mDir/37N-FLAG-IP-mouse_STARAligned.sortedByCoord.out.bam \
              -o $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.bam \
              1:1 2:1 3:1 4:1 5:1 6:1 7:1 8:1 9:1 10:1 11:1 12:1 13:1 14:1 15:1 16:1 17:1 18:1 19:1 X:1 Y:1
samtools index $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.bam

samtools index $mDir/Mock-Total-DI-mouse_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $mDir/Mock-Total-DI-mouse_STARAligned.sortedByCoord.out.bam \
              -o $mDir/Mock-Total-DI-mouse_STARAligned.filtered.bam \
              1:1 2:1 3:1 4:1 5:1 6:1 7:1 8:1 9:1 10:1 11:1 12:1 13:1 14:1 15:1 16:1 17:1 18:1 19:1 X:1 Y:1
samtools index $mDir/Mock-Total-DI-mouse_STARAligned.filtered.bam

# vaccinia
samtools index $vDir/Mock-FLAG-IP-vaccinia_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $vDir/Mock-FLAG-IP-vaccinia_STARAligned.sortedByCoord.out.bam \
              -o $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.bam \
              NC_006998.1:1
samtools index $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.bam

samtools index $vDir/Mock-Total-DI-vaccinia_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $vDir/Mock-Total-DI-vaccinia_STARAligned.sortedByCoord.out.bam \
              -o $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.bam \
              NC_006998.1:1
samtools index $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.bam

samtools index $vDir/37N-FLAG-IP-vaccinia_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $vDir/37N-FLAG-IP-vaccinia_STARAligned.sortedByCoord.out.bam \
              -o $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.bam \
              NC_006998.1:1
samtools index $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.bam

samtools index $vDir/37N-Total-DI-vaccinia_STARAligned.sortedByCoord.out.bam
samtools view -hb -f 2 $vDir/37N-Total-DI-vaccinia_STARAligned.sortedByCoord.out.bam \
              -o $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.bam \
              NC_006998.1:1
samtools index $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.bam

```

#### remove pcr duplicates

Removing duplicates can be done with UMI if UMI barcodes were used, but for this project we're able to just use picard. The exact tool to use is MarkDuplicates, so it's necessary to specify the option REMOVE_DUPLICATES=TRUE.

```
module load picard-2.26.2-gcc-12.1.0

# mouse

picard MarkDuplicates \
      I=$mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.bam \
      O=$mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.dupRm.bam \
      M=marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$mDir/Mock-Total-DI-mouse_STARAligned.filtered.bam \
      O=$mDir/Mock-Total-DI-mouse_STARAligned.filtered.dupRm.bam \
      M=Mock-Total-DI-mouse.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$mDir/37N-FLAG-IP-mouse_STARAligned.filtered.bam \
      O=$mDir/37N-FLAG-IP-mouse_STARAligned.filtered.dupRm.bam \
      M=37N-FLAG-IP-mouse.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$mDir/37N-Total-DI-mouse_STARAligned.filtered.bam \
      O=$mDir/37N-Total-DI-mouse_STARAligned.filtered.dupRm.bam \
      M=37N-Total-DI-mouse.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

# vaccinia

picard MarkDuplicates \
      I=$vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.bam \
      O=$vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.bam \
      M=$vDir/Mock-FLAG-IP-vaccinia.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.bam \
      O=$vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.dupRm.bam \
      M=$vDir/Mock-Total-DI-vaccinia.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.bam \
      O=$vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.bam \
      M=$vDir/37N-FLAG-IP-vaccinia.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE

picard MarkDuplicates \
      I=$vDir/37N-Total-DI-vaccinia_STARAligned.filtered.bam \
      O=$vDir/37N-Total-DI-vaccinia_STARAligned.filtered.dupRm.bam \
      M=$vDir/37N-Total-DI-vaccinia.marked_dup_metrics.txt \
      REMOVE_DUPLICATES=TRUE
```

#### r1 retrieval

Here, some of the code would be different if we were using eCLIP data (where we'd keep R2 instead of R1), and if we had individual replicates that we would need to pool before continuing. However, we're just going to retrieve R1 reads here.

```
# mouse

samtools view -hb -f 130 $mDir/37N-Total-DI-mouse_STARAligned.filtered.dupRm.bam \
              -o $mDir/37N-Total-DI-mouse_STARAligned.filtered.dupRm.R1.bam
samtools index $mDir/37N-Total-DI-mouse_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.dupRm.bam \
              -o $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.dupRm.R1.bam
samtools index $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $mDir/Mock-Total-DI-mouse_STARAligned.filtered.dupRm.bam \
              -o $mDir/Mock-Total-DI-mouse_STARAligned.filtered.dupRm.R1.bam
samtools index $mDir/Mock-Total-DI-mouse_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.dupRm.bam \
              -o $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.dupRm.R1.bam
samtools index $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.dupRm.R1.bam

# vaccinia
samtools view -hb -f 130 $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.dupRm.bam \
              -o $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.dupRm.R1.bam
samtools index $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.bam \
              -o $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.R1.bam
samtools index $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.dupRm.bam \
              -o $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.dupRm.R1.bam
samtools index $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.dupRm.R1.bam

samtools view -hb -f 130 $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.bam \
              -o $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.R1.bam
samtools index $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.dupRm.R1.bam
```

#### bam qc

Finally, we can run fastqc on each sample to check the quality of the data. I'm expecting to have far fewer reads than initially, and obviously the mock vaccinia samples will be empty. But overall I hope the quality will be good!

`fastqc $mDir/*.dupRm.R1.bam`
`fastqc $vDir/*.dupRm.R1.bam`

*QC Note*
Unfortunately, it looks like the picard MarkDuplicates tool, even with the Remove Duplicates option selected, didn't remove all the duplicates - the fastq files here are showing very high percentages of duplicate reads. Because of this, I'm going to try using samtool rmdup instead. And... rmdup is apparently deprecated, so I will try samtools markdup instead, and use the option to remove duplicate reads.

```
# vaccinia
samtools collate -o $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.bam
samtools fixmate -m $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.fixmate.bam
samtools sort -o $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                 $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.fixmate.bam
samtools markdup -r -f $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
                    $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                    $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.bam

samtools collate -o $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.bam
samtools fixmate -m $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.bam
samtools sort -o $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                 $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.bam
samtools markdup -r -f $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
                    $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                    $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.bam

samtools collate -o $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.bam
samtools fixmate -m $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.fixmate.bam
samtools sort -o $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                 $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.fixmate.bam

# killed, needs more than 16G memory
samtools markdup -r -f $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                    $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.bam

samtools collate -o $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.bam
samtools fixmate -m $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.namecollate.bam \
                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.bam
samtools sort -o $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                 $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.bam

# killed, needs more than 16G memory
samtools markdup -r -f $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.stats.txt \
                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.fixmate.coordSort.bam \
                    $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.bam

# mouse
samtools collate -o $mDir/Mock-Total-DI-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.bam
samtools fixmate -m $mDir/Mock-Total-DI-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.fixmate.bam
samtools sort -o $mDir/Mock-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                 $mDir/Mock-Total-DI-mouse_STARAligned.filtered.fixmate.bam
samtools markdup -r -f $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.stats.txt \
                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                    $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.bam

samtools collate -o $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.bam
samtools fixmate -m $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.fixmate.bam
samtools sort -o $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                 $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.fixmate.bam
samtools markdup -r -f $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.stats.txt \
                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                    $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam

samtools collate -o $mDir/37N-Total-DI-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.bam
samtools fixmate -m $mDir/37N-Total-DI-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.fixmate.bam
samtools sort -o $mDir/37N-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                 $mDir/37N-Total-DI-mouse_STARAligned.filtered.fixmate.bam

samtools markdup -r -f $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.stats.txt \
                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                    $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.bam

samtools collate -o $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.bam
samtools fixmate -m $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.namecollate.bam \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.fixmate.bam
samtools sort -o $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                 $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.fixmate.bam
samtools markdup -r -f $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.stats.txt \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.fixmate.coordSort.bam \
                    $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam

```

Then the R1 retrieval needs to be redone:

```
# mouse

samtools view -hb -f 130 $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.bam \
              -o $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam
samtools index $mDir/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam \
              -o $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam
samtools index $mDir/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.bam \
              -o $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam
samtools index $mDir/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.bam \
              -o $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam
samtools index $mDir/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam

# vaccinia
samtools view -hb -f 130 $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.bam \
              -o $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam
samtools index $vDir/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.bam \
              -o $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam
samtools index $vDir/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.bam \
              -o $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam
samtools index $vDir/Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam

samtools view -hb -f 130 $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.bam \
              -o $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam
samtools index $vDir/Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam
```

## pureclip

### basic mode

The parameter -iv allows you to select a subset of chromosomes to 'learn the parameters of PureCLIPs HMM'. Learning on a subset of chromosomes can reduce memory consumption and runtime and shouldn't impair the results; however, I probably will try without this option first. The -nt parameter is the number of threads, for parallelization.

```
# example code from pureclip docs

pureclip \
  -i aligned.prepro.R2.bam \
  -bai aligned.prepro.R2.bam.bai \
  -g ref.hg19.fa \
  -iv 'chr1;chr2;chr3;' -nt 10 \
  -o PureCLIP.crosslink_sites.bed
```

### pureclip with control data

Because we have control data - the Total samples - this is a better way for us to analyze the data. This requires the preprocessed bam files from the inputs and IPs to be submitted together.

The `-o` parameter provides a point-by-point binding site bed file, while the `-or` parameter merges binding sites that are within 8bp of each other (this length can be manually adjusted, but I stuck with the default for now).

```
# example code from pureclip docs

pureclip \
  -i aligned.prepro.R2.bam \
  -bai aligned.prepro.R2.bam.bai \
  -g ref.hg19.fa \
  -iv 'chr1;chr2;chr3;' -nt 10 \
  -o PureCLIP.crosslink_sites.cov_inputSignal.bed \
  -ibam input.aligned.prepro.R2.bam \
  -ibai input.aligned.prepro.R2.bam.bai

# code for mouse

pureclip \
  -i Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -bai Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai \
  -g refDir=/data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  -nt 256 \
  -o Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.bed \
  -or Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.bed \
  -ibam Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -ibai Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai

pureclip \
  -i 37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -bai 37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai \
  -g refDir=/data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.dna.primary_assembly.fa \
  -nt 256 \
  -o 37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.bed \
  -or 37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.bed \
  -ibam 37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam \
  -ibai 37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai

# code for vaccinia

pureclip \
  -i Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  -bai Mock-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai \
  -g refDir=refDir=/data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.fna \
  -nt 256 \
  -o Mock.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.bed \
  -or Mock.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.bed \
  -ibam Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  -ibai Mock-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai

pureclip \
  -i 37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  -bai 37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai \
  -g refDir=refDir=/data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.fna \
  -nt 256 \
  -o 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.bed \
  -or 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.bed \
  -ibam 37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  -ibai 37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai
```

PureCLIP is really struggling to complete the analysis for the Mock samples, which I think is because those samples have almost twice as much data. I'm going to take 50% of the BAM file for both the IP and Input and try running it again. Hopefully this will work!

```
module load samtools-1.21-gcc-12.1.0
samtools view -s 0.05 -b -o Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.small.bam Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam
samtools view -s 0.05 -b -o Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.small.bam Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam

samtools index Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.small.bam
samtools index Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.small.bam
```

Unfortunately, running this command with the subsampled data still failed. I'm not sure if there is something going on with the data, since the size shouldn't be an issue at this point unless there just aren't any binding sites in the mock.

I'm going to try running a very small subsample (0.05 percent) as well as a run just using ip without incorporating the input as a control, to see if either of those work. When I used 50% of the BAM files for IP only, excluding the control data, I was able to call peaks, and when I ran the very small subsample including the Input control I was able to call peaks. I'm going to compare how those peaks look in IGV to see if similar results were identified.

### bedtools manipulation

The output bed files from pureclip can be manipulated with bedtools to get human readable and IGV-viewable files.

```
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/pureclip-env

mdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/mouse"
vdir="/data/gencore/analysis_projects/8406252_Bonner/fastq/vaccinia"

#sort -k1,1 -k2,2n \
#     "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

#sort -k1,1 -k2,2n \
#     "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

#sort -k1,1 -k2,2n \
#     "$mdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.run2.merged.bed \
#     > "$mdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed

sort -k1,1 -k2,2n \
     "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.bed \
     > "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.bed


sort -k1,1 -k2,2n \
     "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.small.run.merged.bed \
     > "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.small.run.merged.sorted.bed
```

To identify the genes that are present in the bed file, I'm going to try using bedtools intersect instead. (`bedtools intersect -wa -a file.gtf -b second_file.bed`).

The intersections, after removing the chromosomal inclusion lines, are quite useful. However, at least in the Mock mouse samples where input controls weren't included in the analysis, some of the called peaks correspond to locations where the input RNA also had high alignment. To correct for this, I want to obtain a list of genes where alignment coverage is at least 10 times higher in the IP than in the Input. For this, I'm using `samtools depth` to obtain the alignment coverage values in each sample's BAM file for the positions in the called peak BED file.

I normalized by the total aligned read depth in the BAM files obtained for IP and Input so that accurate depth comparisons could be made, and used that normalized and filtered BED file to intersect with the reference GTF instead of using the merged BED file directly from PureCLIP.

The following is the code for stringently filtering the PureCLIP peak calls by normalized read depth, to obtain a list of unique genes within which a peak was called. I'm going to use the same stringent filter for both the mock samples (called without input control data) and the 37N samples (called with input control data); if the peak calls with the input control data are more accurate, the only difference should be that a smaller percentage of them are eliminated here.

For the Mock data, there were 944 regions in the BED file and 40 unique genes following normalization and filtering.

```
# calculate depth
samtools depth -b "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.bed \
  -a "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.bam \
  > "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.txt

# identify scaling constant
samtools view -c "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.bam # 899215

# normalize depths
awk '{$4 = $3 / 0.899215} 1' "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.txt \
  | column -t > "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt

# calculate depth
samtools depth -b "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.bed \
  -a "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.bam \
  > "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.txt

# identify scaling constant
samtools view -c "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.bam # 813448

# normalize depths
awk '{$4 = $3 / 0.813448} 1' "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.txt \
  | column -t > "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt

# filter normalized files by relative depths
paste "$mdir"/Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt  \
      "$mdir"/Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt \
      > "$mdir"/Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt

# create a bed file containing only the regions where the normalized depth for the IP samples
# is at least 10 times the normalized depth for the Input samples

awk '$8 >= $4 * 10 {print $1,$2}' "$mdir"/Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt \
  | awk '$3 = $2+1' | tr ' ' \\t > "$mdir"/Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.bed

# rerun the intersection using the BED file filtered by relative normalized read depths
bedtools intersect \
  -wa -b "$mdir"/Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.bed \
  -a /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.114.sorted.gff3 \
  > "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.gtf

# remove the chromosome regions from the gtf file
# keep only the unique regions in the gtf file
grep -v 'ID=chromosome' "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.gtf \
  | sort | uniq > "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.rmChr.gtf

# identify only the unique genes involved in identified peaks
grep 'ID=gene' "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.rmChr.gtf \
  > "$mdir"/Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf
```

For the 37N mouse data there were 209 regions in the original BED file from PureCLIP, and only 6 genes remained when filtered in this way (from 34 total annotations). Interestingly, changing the cutoff value from IP at least 10X Input reads to IP only 2X Input reads didn't actually change the number of unique annotations or genes.

```
# calculate depth
samtools depth -b "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam \
  > "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.depth.txt

# identify scaling constant
samtools view -c "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam # 695556

# normalize depths
awk '{$4 = $3 / 0.695556} 1' "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.depth.txt \
  | column -t > "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# calculate depth
samtools depth -b "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam \
  > "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.depth.txt

# identify scaling constant
samtools view -c "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam # 863894

# normalize depths
awk '{$4 = $3 / 0.863894} 1' "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.depth.txt \
  | column -t > "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# filter normalized files by relative depths
paste "$mdir"/37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt  \
      "$mdir"/37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt \
      > "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# create a bed file containing only the regions where the normalized depth for the IP samples
# is at least 10 times the normalized depth for the Input samples

awk '$8 >= $4 * 10 {print $1,$2}' "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt \
  | awk '$3 = $2+1' | tr ' ' \\t > "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.bed

# rerun the intersection using the BED file filtered by relative normalized read depths
bedtools intersect \
  -wa -b "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.bed \
  -a /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.114.sorted.gff3 \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.gtf

# remove the chromosome regions from the gtf file
# keep only the unique regions in the gtf file
grep -v 'ID=chromosome' "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.gtf \
  | sort | uniq > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf

# identify only the unique genes involved in identified peaks
grep 'ID=gene' "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf


# create a bed file containing only the regions where the normalized depth for the IP samples
# is at least twice the normalized depth for the Input samples

awk '$8 >= $4 * 2 {print $1,$2}' "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt \
  | awk '$3 = $2+1' | tr ' ' \\t > "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.2X.bed

# rerun the intersection using the BED file filtered by relative normalized read depths
bedtools intersect \
  -wa -b "$mdir"/37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.2X.bed \
  -a /data/gencore/databases/reference_genomes/mouse/GRCm39/Mus_musculus.GRCm39.114.sorted.gff3 \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.2X.normFiltered.intersect.gtf

# remove the chromosome regions from the gtf file
# keep only the unique regions in the gtf file
grep -v 'ID=chromosome' "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.2X.normFiltered.intersect.gtf \
  | sort | uniq > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.2X.normFiltered.intersect.rmChr.gtf

# identify only the unique genes involved in identified peaks
grep 'ID=gene' "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.2X.normFiltered.intersect.rmChr.gtf \
  > "$mdir"/37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.2X.normFiltered.intersect.rmChr.onlyGenes.gtf
```

I also ran this normalization and filtering for the 37N vaccinia data. Obviously the Mock vaccinia alignment data is non-existent, so no comparisons can be made after identification of peaks.

```
# calculate depth
samtools depth -b "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  > "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.depth.txt

# identify scaling constant
samtools view -c "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam # 971024

# normalize depths
awk '{$4 = $3 / 0.971024} 1' "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.depth.txt \
  | column -t > "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# calculate depth
samtools depth -b "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam \
  > "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.depth.txt

# identify scaling constant
samtools view -c "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam # 1839647

# normalize depths
awk '{$4 = $3 / 1.839647} 1' "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.depth.txt \
  | column -t > "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# filter normalized files by relative depths
paste "$vdir"/37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt  \
      "$vdir"/37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt \
      > "$vdir"/37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt

# create a bed file containing only the regions where the normalized depth for the IP samples
# is at least 10 times the normalized depth for the Input samples

awk '$8 >= $4 * 10 {print $1,$2}' "$vdir"/37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt \
  | awk '$3 = $2+1' | tr ' ' \\t > "$vdir"/37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.bed

# rerun the intersection using the BED file filtered by relative normalized read depths
bedtools intersect \
  -wa -b "$vdir"/37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.bed \
  -a /data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.gtf \
  > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.gtf

# remove the chromosome regions from the gtf file
# keep only the unique regions in the gtf file
grep -v 'ID=chromosome' "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.gtf \
  | sort | uniq > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf

# identify only the unique genes involved in identified peaks
grep 'transcript_id ""' "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf \
  > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf
```

Twenty genes were identified after filtering, but I was curious to know how many had been predicted by PureCLIP originally. The following code showed that 138 genes had originally been identified as potential binding sites, which supports possible non-specific binding of the protein of interest to vaccinia RNAs.

```
bedtools intersect \
  -wa -b "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.bed \
  -a /data/gencore/databases/reference_genomes/vaccinia/GCF_000860085.1_ViralProj15241_genomic.gtf \
  > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.intersect.gtf

grep 'transcript_id ""' "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.intersect.gtf \
  | sort | uniq > "$vdir"/37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.intersect.onlyGenes.gtf
```

## Rename Final Output Files

The file names generated in this workflow are very long and unwieldy, so I'm going to copy the output files for James into renamed files (so I don't lose the original naming for my backup data) and keep track of the renaming here.

```
cp 37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam 37N-FLAG-IP-vaccinia.processed.bam
cp 37N-FLAG-IP-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai 37N-FLAG-IP-vaccinia.processed.bam.bai
cp 37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam 37N-Total-DI-vaccinia.processed.bam
cp 37N-Total-DI-vaccinia_STARAligned.filtered.samDupRm.R1.bam.bai 37N-Total-DI-vaccinia.processed.bam.bai
cp 37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.bed 37N.vaccinia.depthNormalized.bed
cp 37N.vaccinia_STARAligned.filtered.samDupRm.R1.depth.normalized.txt 37N.vaccinia.depthNormalized.txt
cp 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf 37N.vaccinia.crosslinks.gtf
cp 37N.vaccinia.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf 37N.vaccinia.crosslinks.genes.gtf

cp 37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam 37N-FLAG-IP-mouse.processed.bam
cp 37N-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai 37N-FLAG-IP-mouse.processed.bam.bai
cp 37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam 37N-Total-DI-mouse.processed.bam
cp 37N-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai 37N-Total-DI-mouse.processed.bam.bai
cp 37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.bed 37N.mouse.depthNormalized.bed
cp 37N.mouse_STARAligned.filtered.samDupRm.R1.depth.normalized.txt 37N.mouse.depthNormalized.txt
cp 37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.gtf 37N.mouse.crosslinks.gtf
cp 37N.mouse.PureCLIP.crosslink_sites.cov_inputSignal.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf 37N.mouse.crosslinks.genes.gtf

cp Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam Mock-FLAG-IP-mouse.processed.bam
cp Mock-FLAG-IP-mouse_STARAligned.filtered.samDupRm.R1.bam.bai Mock-FLAG-IP-mouse.processed.bam.bai
cp Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam Mock-Total-DI-mouse.processed.bam
cp Mock-Total-DI-mouse_STARAligned.filtered.samDupRm.R1.bam.bai Mock-Total-DI-mouse.processed.bam.bai
cp Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.bed Mock.mouse.depthNormalized.bed
cp Mock.mouse_STARAligned.filtered.samDupRm.R1.subsampled.depth.normalized.txt Mock.mouse.depthNormalized.txt
cp Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.rmChr.gtf Mock.mouse.crosslinks.gtf
cp Mock.mouse.PureCLIP.crosslink_sites.cov_inputSignal.subsampled.runNC.merged.sorted.normFiltered.intersect.rmChr.onlyGenes.gtf Mock.mouse.crosslinks.genes.gtf
```
