# Code Used for Differential Gene Expression Analysis for iLab 8424037

## Module A: Alignment and Quantification

These samples are bovine, and since I don't have STAR indexes made I'll be running genome generate.

```
refDir="/data/gencore/databases/reference_genomes/bovine/bos-taurus-UMD3.1"
fasta="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.fna"
gff="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gff"
gtf="$refDir"/"GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gtf"

cd "$refDir"

gffread "$gff" -T -o "$gtf"

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

While the indexes were building, I removed ribosomal RNA using bbduk, then used cutadapt and trimmomatic to remove any remaining adapter sequences or very low quality reads. The sequence I used for the ribosomal RNA can be found at https://www.ncbi.nlm.nih.gov/nuccore/DQ222453.1?report=fasta, and includes the *Bos taurus* 18S, ITS1, 5.8S, ITS2, and 28S sequences.

```
module load mamba/latest

scriptsDir="/data/gencore/shared_scripts/RNAseq/old-scripts"
refSeq="/data/gencore/databases/reference_genomes/bovine/bos.taurus.rRNA.fna"

fastqDir="/data/gencore/analysis_projects/8424037_Molehin/fastq"
bbstatsDir="/data/gencore/analysis_projects/8424037_Molehin/bbstats"

mkdir -p "$bbstatsDir"

source activate /data/biocore/programs/conda-envs/bb-env

cd $fastqDir

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  echo $i
  bbduk.sh in="$i"_SRN_L001_R1_001.fastq.gz in2="$i"_SRN_L001_R2_001.fastq.gz \
           outm1="$bbstatsDir"/ribo-"$i"_SRN_L001_R1_001.fastq outm2="$bbstatsDir"/ribo-"$i"_SRN_L001_R2_001.fastq \
           outu1="$bbstatsDir"/nonribo-"$i"_SRN_L001_R1_001.fastq outu2="$bbstatsDir"/nonribo-"$i"_SRN_L001_R2_001.fastq \
           ref=$refSeq stats="$bbstatsDir"/"$i"_bbdukstats.txt
done;

source deactivate

cd $bbstatsDir
chmod -R g+w *
```

The remaining non-ribosomal reads were then filtered and trimmed.

```
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/cutadapt/

mkdir -p /data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/cutadapt
cd /data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo

adapters="/data/gencore/databases/trimmomatic/PolyAndIllumina.fa"

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  cutadapt -a file:"$adapters" -A file:"$adapters" \
         -m 30 -q 20 \
         -o "$i"-cut_SCT_L001_R1_001.fastq.gz \
         -p "$i"-cut_SCT_L001_R2_001.fastq.gz \
         "$i"_S*_L001_R1_001.fastq.gz "$i"_S*_L001_R2_001.fastq.gz
done

mv *cut* ./cutadapt

cd ./cutadapt

module load trimmomatic-0.39-fs

for i in $(find ./ -type f -name "*.fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1)
do
  trimmomatic PE "$i"_SCT_L001_R1_001.fastq.gz "$i"_SCT_L001_R2_001.fastq.gz \
            "$i"_SQP_L001_R1_001.fastq.gz "$i"_SUN_L001_R1_001.fastq.gz \
            "$i"_SQP_L001_R2_001.fastq.gz "$i"_SUN_L001_R2_001.fastq.gz \
            CROP:150 SLIDINGWINDOW:4:15 MINLEN:36
done

mkdir -p ../cut-paired-filter0415
mkdir -p ../cut-unpaired-filter0415

mv *SQP* ../cut-paired-filter0415/
mv *SUN* ../cut-unpaired-filter0415/

cd ../
chmod -R g+w *
```

Because this is mammalian data, I used the default recommended STAR parameters for alignment.

```
sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8424037_Molehin/bbstats/nonribo/cut-paired-filter0415/fastq \
  -i /scratch/kawoodbu/8424037_Molehin_UMD-nonribo_default \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-default.txt \
  -r /data/gencore/databases/reference_genomes/bovine/bos-taurus-UMD3.1 \
  -a /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut \
  -q /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

This worked quite well! Because of the two different aims for this experiment, I separated out the samples to quantify and run differential analysis on, so the modeling formula will only take into consideration samples from the appropriate portion of the experiment.

```
# aim 1
for i in $(find /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM1 -type f -name "*.bam" | \
            while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq )
do
echo $i
bash /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A2.stringtie-quant.sh \
  -s $i -a /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM1 \
  -q /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM1 \
  -g /data/gencore/databases/reference_genomes/bovine/bos-taurus-UMD3.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gtf
done

bash /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A \
  -q /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM1 \
  -a /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM1

# aim2
for i in $(find /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM2 -type f -name "*.bam" | \
            while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq )
do
echo $i
bash /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A2.stringtie-quant.sh \
  -s $i -a /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM2 \
  -q /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM2 \
  -g /data/gencore/databases/reference_genomes/bovine/bos-taurus-UMD3.1/GCF_000003055.6_Bos_taurus_UMD_3.1.1_genomic.gtf
done

bash /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A \
  -q /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM2 \
  -a /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-alignment-default-nonribo-cut/AIM2
```

## Module B: Differential Expression

```
sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-differentials-default-nonribo-cut/AIM1 \
  -g /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM1/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8424037_Molehin/AIM1-comparisons.csv \
  -s "deseq2 edger noiseq"

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-differentials-default-nonribo-cut/AIM2 \
  -g /data/gencore/analysis_projects/8424037_Molehin/bos-UMD-quants-default-nonribo-cut/AIM2/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8424037_Molehin/AIM2-comparisons.csv \
  -s "deseq2 edger noiseq"
```

The next step is to merge the results into a summary output of the three tools, and create some Venn diagrams to visualize the overlap in DEG calls - basically the following command run for both comparison groups. I chose to run merging without a log2 fold change cutoff, but the csv files can be manually filtered if only genes with a higher log2 fold change in expression are of interest.

```
source activate /data/biocore/programs/mamba-envs/biocore-rna
python /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/merge_de_both.py /data/gencore/analysis_projects/8424037_Molehin/AIM1-comparisons.csv 0
python /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/merge_de_both.py /data/gencore/analysis_projects/8424037_Molehin/AIM2-comparisons.csv 0
```

## Module C: Functional Enrichment

For this project, I am using PANTHER to identify functional terms that are over-represented in the lists of DEGs from each comparison, since a gene ontology file is not available for the reference genome, and since *Bos taurus* is available in the PANTHER database.

To make the input files for PANTHER (these are for AIM2 but I did the same type of command for the comparisons in AIM1) I used the following code:

```
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.Calf_vs_Cow.up.0.0.stats.csv > calf_vs_cow_up0_pantherInput.txt
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.Calf_vs_Cow.down.0.0.stats.csv > calf_vs_cow_down0_pantherInput.txt
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.Crypto_vs_UN.up.0.0.stats.csv > crypto_vs_un_up0_pantherInput.txt
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.Crypto_vs_UN.down.0.0.stats.csv > crypto_vs_un_down0_pantherInput.txt
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.C_vs_B.up.0.0.stats.csv > c_vs_b_up0_pantherInput.txt
awk -F "[|,]" -v OFS='\t' '{ print $2 }' merged_deg.C_vs_B.down.0.0.stats.csv > c_vs_b_down0_pantherInput.txt
```

I used the functional overrepresentation test in PANTHER (version 19.0, analysis release 20240807) with Fisher's exact test and the Bonferroni correction for multiple testing. The resulting functional term lists were visualized with dot plots in R. The gene ontology (GO) database used was the version released 2025-03-16.
