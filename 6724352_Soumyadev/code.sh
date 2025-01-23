# create annotation reference from Massilia GTF file
awk -F '\t' '{ print $9 }' Massilia.MET4.gff | awk -F ';' '{ print $1"\t"$2 }' | sed 's/Name=//g' | sed 's/ID=//g' > massilia.gene.info.txt

# cpu-align script (testing, unwrapped)
for i in $(find "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/cpu-align.sh \
    --sampleID "$i" \
    --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
    --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
    --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-default.txt" \
    --refDir "/data/gencore/databases/reference_genomes/massilia" \
    --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment"
done

# some of the samples failed for a lack of memory, so I've increased the memory parameter in the script
sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A1.cpu-align.sh \
  --sampleID "Massila-plus-A" \
  --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
  --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
  --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-default.txt" \
  --refDir "/data/gencore/databases/reference_genomes/massilia" \
  --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment"

sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A1.cpu-align.sh \
  --sampleID "Massila-plus-B" \
  --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
  --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
  --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-default.txt" \
  --refDir "/data/gencore/databases/reference_genomes/massilia" \
  --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment"

sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A1.cpu-align.sh \
  --sampleID "Massila-plus-C" \
  --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
  --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
  --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-default.txt" \
  --refDir "/data/gencore/databases/reference_genomes/massilia" \
  --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/default-alignment"

#### all samples had poor alignment rates with default params, especially the plus samples
# trying again with bac lowqual
for i in $(find "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A1.cpu-align.sh \
    --sampleID "$i" \
    --fastqDir "/data/gencore/analysis_projects/6724352_Soumyadev/fastq" \
    --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev" \
    --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-bacterial-lowqual.txt" \
    --refDir "/data/gencore/databases/reference_genomes/massilia" \
    --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-alignment"
done

# alignment really isn't looking good, with almost 50% of reads aligning to multiple loci
# here, create fasta file for massilia ribsome sequence to check rRNA content with bbduk
touch massilia.rRNA.bed
vim insert: contig_1 2226552 2227693

module load bedtools2-2.30.0-gcc-11.2.0
bedtools getfasta -fi Massilia_assembly.fasta -bed massilia.rRNA.bed -fo massilia.rRNA.fa

# bbduk was fine, and the alignment stats look similar to massilia-plus samples from his previous run
# so I'm just going to continue.

# run stringtie quant once all 6 samples are complete:
for i in $(find "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-alignment" -type f -name "*.final.out" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A2.stringtie-quant.sh \
    --sampleID "$i" \
    --alignmentDir "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-alignment" \
    --quantDir "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-stringtie" \
    --refGTF "/data/gencore/databases/reference_genomes/massilia/Massilia.MET4.gff"
done

# merge stringtie results to create gene count matrix
sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A3.merge-quant.sh \
  --scriptDir "/data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024" \
  --quantDir "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-stringtie" \
  --alignmentDir "/data/gencore/analysis_projects/6724352_Soumyadev/bac-lowqual-alignment"

# run deseq R scripts
sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/B1.DEG.rscripts.sh --script noiseq \
      --directory "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-DEGs" \
      --geneMatrix "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-stringtie/transcript_count_matrix.csv" \
      --comparisons "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/comparisons.csv"

sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/B1.DEG.rscripts.sh --script edger \
      --directory "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-DEGs" \
      --geneMatrix "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-stringtie/transcript_count_matrix.csv" \
      --comparisons "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/comparisons.csv"

sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/B1.DEG.rscripts.sh --script deseq2 \
      --directory "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-DEGs" \
      --geneMatrix "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/bac-lowqual-stringtie/transcript_count_matrix.csv" \
      --comparisons "/data/gencore/analysis_projects/completed_projects/6724352_Soumyadev/comparisons.csv"

# run this script in directory containing tabular output from the three DEG tools to merge at different fold changes (0, 1, 2)
python /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/merge_de_both.py ../../comparisons.csv 2.0

# annotate the tables with the data in the gene info text file in the reference genome directory
bash /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/C2.annotate-DEGs.sh \
  --infoFile /data/gencore/databases/reference_genomes/massilia/massilia.gene.info.txt \
  --inputDirectory /Volumes/Gencore/sftp/f_garciapichel/6724352_Soumyadev_RNA/3.differential-expression-corrected/edgeR-deglists

#### CODE WITH TRIMMED READS ####
# wondering if trimming adapter content would change things significantly
sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/0.trimmomatic.sh \
  --directory "/data/gencore/analysis_projects/5663018_Masmudur_RNA/project_6/fastq/fastq" \
  --outputDir "/data/gencore/analysis_projects/5663018_Masmudur_RNA/project_6/trimmed-fastq"

# trimming changed the read counts by a lot, so lets continue
for i in $(find "/data/gencore/sftp/f_garciapichel/6724352_Soumyadev_RNA/adapter-trimmed-analysis/1.sequencing/fastq" -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d '_' -f 1 | sort | uniq)
do
  echo "$i"
  sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/A1.cpu-align.sh \
    --sampleID "$i" \
    --fastqDir "/data/gencore/sftp/f_garciapichel/6724352_Soumyadev_RNA/adapter-trimmed-analysis/1.sequencing/fastq" \
    --intermediateDir "/scratch/kawoodbu/6724352_Soumyadev_trimmed" \
    --starParams "/data/gencore/analysis_projects/6724352_Soumyadev/star-params-bacterial-lowqual.txt" \
    --refDir "/data/gencore/databases/reference_genomes/massilia" \
    --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/trimmed-bac-lowqual-alignment"
done

bash /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/2.alignment-wrapper.sh \
        -f /data/gencore/sftp/f_garciapichel/6724352_Soumyadev_RNA/adapter-trimmed-analysis/1.sequencing/fastq \
        -i /scratch/kawoodbu/6724352_Soumyadev_trimmed \
        -p /data/gencore/analysis_projects/6724352_Soumyadev/star-params-bacterial-lowqual.txt \
        -r /data/gencore/databases/reference_genomes/massilia \
        -o /data/gencore/analysis_projects/6724352_Soumyadev/trimmed-bac-lowqual-alignment \
        -q /data/gencore/analysis_projects/6724352_Soumyadev/trimmed-bac-lowqual-stringtie \
        -s /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024

### based on some forums I think I can just trim adapters and not QC, and try to preserve more reads.
## still keeping the min length for trimming at 100 though
## this loses a LOT of reads, but I'm still going to run alignment to see how the stats compare
## if only unmapped reads are lost, it'll be worth it.
## Nope! Not worth it to trim, for either this project or Masmudur's, so probably not in general.
sbatch /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/0.trimmomatic.sh \
  --directory "/data/gencore/sftp/f_garciapichel/6724352_Soumyadev_RNA/untrimmed-analysis/1.sequencing/fastq" \
  --outputDir "/data/gencore/analysis_projects/6724352_Soumyadev/trimmed-fastq"

bash /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024/2.alignment-wrapper.sh \
        -f //data/gencore/analysis_projects/6724352_Soumyadev/trimmed-fastq/adapter-trimmed/fastq \
        -i /scratch/kawoodbu/6724352_Soumyadev_trimmed \
        -p /data/gencore/analysis_projects/6724352_Soumyadev/star-params-bacterial-lowqual.txt \
        -r /data/gencore/databases/reference_genomes/massilia \
        -o /data/gencore/analysis_projects/6724352_Soumyadev/adpt-trim-bac-lowqual-alignment \
        -q /data/gencore/analysis_projects/6724352_Soumyadev/adpt-trim-bac-lowqual-stringtie \
        -s /data/gencore/shared_scripts/RNAseq/updated-pathway-jan2024
