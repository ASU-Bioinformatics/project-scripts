Great, I have a typo in the environment name lol.

```
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/ioseq-env
```

### Map Reads to Reference Genome

Although I mapped first, using the original bam files from the Revio, I think this is still what I needed to use as input. The `collapse` tool will then integrate this mapped bam with the refined flnc bam. (I'm not sure where the clustered bam shows up again).

Actually, it looks like I can use the clustered bam file as input for the alignment - maybe as an alternative? I guess I can try both.

```
pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc03_5p--IsoSeqX_3p.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc03_5p--IsoSeqX_3p.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc04_5p--IsoSeqX_3p.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc04_5p--IsoSeqX_3p.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc11_5p--IsoSeqX_3p.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc11_5p--IsoSeqX_3p.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc12_5p--IsoSeqX_3p.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/umaydis-aligned.IsoSeqX_bc12_5p--IsoSeqX_3p.bam
```

### Refine...
So, I did the read mapping prematurely. To clean up the data, I'll really need to start with refining the file (remove polyA tails and artificial concatemers). I'm trying the first one with the --require-polya parameter because I am not sure whether these samples have poly(A) tails or not...

Running with the polyA parameter detects a polyA tail on the majority of reads, so I'm going to use that parameter so those tails can be trimmed.

```
isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc03_5p--IsoSeqX_3p.bam \
              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p.polya.flnc.bam

isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc04_5p--IsoSeqX_3p.bam \
              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p.polya.flnc.bam

isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc11_5p--IsoSeqX_3p.bam \
              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p.polya.flnc.bam

isoseq refine /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/m84082_241014_083242_s2.hifi_reads.IsoSeqX_bc12_5p--IsoSeqX_3p.bam \
              /data/gencore/analysis_projects/8961526_Dunn/primers.fasta --require-polya \
              /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p.polya.flnc.bam
```

### Cluster Isoforms

```
isoseq cluster2 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.polya.flnc.bam \
                /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.clustered.bam

isoseq cluster2 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.polya.flnc.bam \
                /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.clustered.bam

isoseq cluster2 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.polya.flnc.bam \
                /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.clustered.bam

isoseq cluster2 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.polya.flnc.bam \
                /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.clustered.bam                
```

### Map Clustered Reads to Reference

```
pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.clustered.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.aligned.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.clustered.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.aligned.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.clustered.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.aligned.bam

pbmm2 align --preset ISOSEQ \
  --sort /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.clustered.bam \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.aligned.bam

```

### Collapse Into Unique Isoforms

This step identifies all the unique isoforms and provides counts for each one.

```
isoseq collapse --do-not-collapse-extra-5exons \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.aligned.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.polya.flnc.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.gff

 isoseq collapse --do-not-collapse-extra-5exons \
   /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.aligned.bam \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.polya.flnc.bam \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.gff

isoseq collapse --do-not-collapse-extra-5exons \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.aligned.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.polya.flnc.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.gff

isoseq collapse --do-not-collapse-extra-5exons \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.aligned.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.polya.flnc.bam \
 /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.gff
```

### Pigeon Prepare

Both the reference files and transcript GFFs from `isoseq collapse` need to be prepared. I think I can just list them all in one command:

```
pigeon prepare \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.gff \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.gff \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.gff \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.gff \

pigeon prepare /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.gtf
```

Unfortunately, the reference GTF/GFF file is not correctly formatted for *pigeon prepare*, so I'm going to have to figure out a way to adapt the file to meet their requirements. Most of the formatting was already correct, but the annotations didn't contain the `gene_name` field. It took a lot of trial and error, but to get a compatible file I added a `gene_name` field that was identical to the `gene_id` field, retained only the `gene`, `transcript`, and `exon` entries, removed all the `transcript_id` fields from the `gene` entries because they were empty, and made sure all the metadata header lines began with two hashtags instead of one.

### Pigeon Classify, Filter, and Report

This step classifies the isoforms into categories (full splice match, incomplete splice match, novel in catalog, novel out of catalog, antisense, intron, genomic, and intergenic).

```
pigeon classify \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.sorted.gff \
  /data/gencore/analysis_projects/8961526_Dunn/umaydis.edited.sorted.gtf \
  /data/gencore/databases/reference_genomes/umaydis/GCF_000328475.2/GCF_000328475.2_Umaydis521_2.0_genomic.fna \
  --fl /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.flnc_count.txt
```

Then we filter out artefacts (potentially poor quality data - primarily anything non-canonical) and report the saturation for it.

```
pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt

pigeon filter /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.txt --isoforms /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.sorted.gff

pigeon report --exclude-singletons /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
  /data/gencore/analysis_projects/8961526_Dunn/IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_saturation.txt
```

## Custom Work

### Identify Spliced Transcripts

I'm not really sure of the best way to do this, but I'm starting by filtering out all of the mono-exonic transcripts, since there isn't any splicing going on there. However, some of those may be variants of a spliced transcript so I'm thinking I may need to go about it a different way. Potentially, I can examine any other isoforms belonging to the same gene as the non-mono-exonic isoforms? There are also multiple mono-exonic transcripts for many of the genes...

Ok I'm working on a format for comparing isoforms between samples. To test it, I'm making versions of all the necessary input files that just go through the first three genes in the classification text files (UMAG_06480, UMAG_10945, and UMAG_10944). Conveniently, this corresponds to PB.1-3 isoform IDs. Inconveniently, it covers very few unique junctions.

Constants:
sample id 1: bc03   total FL reads: 1785168
sample id 1: bc04   total FL reads: 3997081
sample id 1: bc11   total FL reads: 3907798
sample id 1: bc12   total FL reads: 4559134

```
# subset classification files and grab relevant columns to decrease memory load on Python later
head -n 15 IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
| awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' > comparisons/bc03_classification.txt

head -n 32 IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
| awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' > comparisons/bc04_classification.txt

head -n 14 IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
| awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' > comparisons/bc11_classification.txt

head -n 15 IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt \
| awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' > comparisons/bc12_classification.txt

# subset junctions files and grab relevant columns to decrease memory load on Python later
head -n 2 IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt \
| awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' > comparisons/bc03_junctions.txt

head -n 2 IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt \
| awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' > comparisons/bc04_junctions.txt

head -n 1 IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt \
| awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' > comparisons/bc11_junctions.txt

head -n 1 IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt \
| awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' > comparisons/bc12_junctions.txt

# subset classification files and grab relevant columns to decrease memory load on Python later
head -n 23 IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.abundance.txt \
| awk -v OFS='\t' '{print $1, $3}' > comparisons/bc03_abundance.txt

head -n 40 IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.abundance.txt \
| awk -v OFS='\t' '{print $1, $3}' > comparisons/bc04_abundance.txt

head -n 21 IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.abundance.txt \
| awk -v OFS='\t' '{print $1, $3}' > comparisons/bc11_abundance.txt

head -n 21 IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.abundance.txt \
| awk -v OFS='\t' '{print $1, $3}' > comparisons/bc12_abundance.txt
```

For simplicity of naming, I just reused the same code as above without subsetting the number of lines to create the full isoform comparison file.

```
# subset classification files and grab relevant columns to decrease memory load on Python later
awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt > comparisons-complete-files/bc03_classification.txt

awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt > comparisons-complete-files/bc04_classification.txt

awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt > comparisons-complete-files/bc11_classification.txt

awk -v OFS='\t' '{print $1, $7, $8, $2, $3, $9, $10, $6, $15, $4, $5}' IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_classification.txt > comparisons-complete-files/bc12_classification.txt

# subset junctions files and grab relevant columns to decrease memory load on Python later
awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' IsoSeqX_bc03_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt > comparisons-complete-files/bc03_junctions.txt

awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' IsoSeqX_bc04_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt > comparisons-complete-files/bc04_junctions.txt

awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' IsoSeqX_bc11_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt > comparisons-complete-files/bc11_junctions.txt

awk -v OFS='\t' '{print $1, $4, $5, $6, $8, $9, $10}' IsoSeqX_bc12_5p--IsoSeqX_3p_classification.filtered_lite_junctions.txt > comparisons-complete-files/bc12_junctions.txt

# subset classification files and grab relevant columns to decrease memory load on Python later
awk -v OFS='\t' '{print $1, $3}' IsoSeqX_bc03_5p--IsoSeqX_3p/IsoSeqX_bc03_5p--IsoSeqX_3p.collapsed.abundance.txt > comparisons-complete-files/bc03_abundance.txt

awk -v OFS='\t' '{print $1, $3}' IsoSeqX_bc04_5p--IsoSeqX_3p/IsoSeqX_bc04_5p--IsoSeqX_3p.collapsed.abundance.txt > comparisons-complete-files/bc04_abundance.txt

awk -v OFS='\t' '{print $1, $3}' IsoSeqX_bc11_5p--IsoSeqX_3p/IsoSeqX_bc11_5p--IsoSeqX_3p.collapsed.abundance.txt > comparisons-complete-files/bc11_abundance.txt

awk -v OFS='\t' '{print $1, $3}' IsoSeqX_bc12_5p--IsoSeqX_3p/IsoSeqX_bc12_5p--IsoSeqX_3p.collapsed.abundance.txt > comparisons-complete-files/bc12_abundance.txt
```

Once I had these sample files, I made the script `merge.py` to rearrange and combine the different relevant fields from within these files. I hardcoded the input values into the script but I will need to go back and update the script if it ends up being a valuable part of the IsoSeq workflow. Also, it looks like something I did in the script that didn't get tested out with the subset data caused some duplication of lines, so I did a `sort | uniq` on the file in bash. I also used `grep` to extract just the multi-exon isoforms, since mono-exon isoforms aren't going to provide much splicing information...

```
grep 'multi-exon' all-isoforms.txt > multi-exon-isoforms.txt
sort multi-exon-isoforms.txt | uniq > multi-exon-isoforms-removedups.txt

sort all-isoforms.txt | uniq > all-isoforms-removedups.txt

sed '/^\0/d' all-isoforms-removedups.txt > all-isoforms-removedups-knowngene.txt
```

In R, I am going to find the abundance of mono-exonic transcript reads for each gene in each sample, as well as the abundance of multi-exonic transcript reads for each gene in each sample.
