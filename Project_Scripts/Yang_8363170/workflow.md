# Code Used for Differential Gene Expression Analysis for iLab 8363170

## Module A: Alignment and Quantification

These samples are dual culture E. coli and Pseudomonas aeruginosa. I don't have STAR indexes made for either species, so I'll be running genome generate twice.
For consistency with past projects from Jiseon, I'm using the CFT073 E. coli strain instead of the K-12 standard.

The first attempt at building the E. coli indexes failed because there was a space in column two for most entries, which threw off STAR's ability to detect which features were CDS. I used the code `sed -i 's/Protein[[:space:]]Homology/ProteinHomology/g' CFT073.genomic.gtf` to remove that extra space.

```
# escherichia coli

refDir="/data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2"
fasta="$refDir"/"CFT073.genomic.fna"
gtf="$refDir"/"CFT073.genomic.gtf"

cd "$refDir"

STAR \
  --runThreadN 24 \
  --runMode genomeGenerate \
  --genomeDir ./ \
  --genomeFastaFiles "$fasta" \
  --sjdbGTFfile "$gtf" \
  --sjdbOverhang 151 \
  --sjdbGTFfeatureExon CDS \
  --genomeSAindexNbases 10 \
  --limitGenomeGenerateRAM 3000000000000
```

```
# pseudomonas aeruginosa

refDir="/data/gencore/databases/reference_genomes/pseudomonas_aeruginosa"
fasta="$refDir"/"GCF_000006765.1_PAO1.genomic.fna"
gtf="$refDir"/"GCF_000006765.1_PAO1.genomic.gtf"

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
The naming for these samples is non-traditional so I'm renaming them.

```
cd /data/gencore/analysis_projects/8363170_Yang/fastq/
rename '_L005' '_SRN_L001' /data/gencore/analysis_projects/8363170_Yang/fastq/*
```

I'm using the GitHub repository for the RNA scripts this time, now that it's connected to Sol!

```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-default.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-default \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-default \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-default.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-default \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-default \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

Unfortunately, most reads showed up as too short or as multimapped for this culture, so I'm going to run a different parameter setting to try to capture those.


```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_ecoli \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

Again here most of the reads were multimapped. STAR will use the multimapped reads, but I'm concerned based on the percentages that many of these reads are binding multiple times to places on the genome that are very similar between species. For this reason I'm going to require slightly longer binding length but reduce the multimapping max number.

```
# escherichia coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_ecoli \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/ecoli/CFT073_GCF_014262945.1_ASM1426294v1/starTry2 \
  -a /data/gencore/analysis_projects/8363170_Yang/ecoli-alignment-bacterial-short30 \
  -q /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30 \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

```
# pseudomonas aeruginosa
# I forgot to change the output folder name for this one so it overwrote the regular bacterial short parameters
# but the results from that were so bad I don't think it really matters

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A
```

There was an unexplainable error for sample Spx21-Ground-Mid-G2AG in the Pseudomonas alignment (it said biocore-rna wasn't a conda environment, when it had worked for all other samples). So, I have to rerun that sample separately and then restart the merge and quant scripts.

```
sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A1.cpu-align.sh \
  -s Spx21-Ground-Mid-G2AG \
  -f /data/gencore/analysis_projects/8363170_Yang/fastq \
  -i /scratch/kawoodbu/8363170_Yang_pseudo \
  -p /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/supplemental_files/star-params-bacterial-short30-3max.txt \
  -r /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa \
  -o /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short

for i in $(find /data/gencore/analysis_projects/8363170_Yang/fastq -type f -name "*fastq.gz" | while read F; do basename $F; done | cut -d "_" -f 1 | sort | uniq)
do
  echo $i
  sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A2.stringtie-quant.sh \
    -s "$i" \
    -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short30 \
    -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short30 \
    -g /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.genomic.gtf
done

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A/A3.merge-quant.sh \
  -s /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_A \
  -q /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short30 \
  -a /data/gencore/analysis_projects/8363170_Yang/paeruginosa-alignment-bacterial-short30
```

## Module B: Differential Expression

After speaking with Jiseon, I set up the comparisons file with only the comparisons Treated v. Untreated and Flight v. Ground. I did also look at some other comparisons, especially between the two cultures, and didn't observe any batch effect.

```
# e. coli

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8363170_Yang/ecoli-differentials-bacterial-short30 \
  -g /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8363170_Yang/comparisons.csv \
  -s "deseq2 edger noiseq"

# e. coli with IDs only, for functional enrichment

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8363170_Yang/ecoli-differentials-idsonly-bacterial-short30 \
  -g /data/gencore/analysis_projects/8363170_Yang/ecoli-quants-bacterial-short30/gene_count_matrix_idsonly.csv \
  -c /data/gencore/analysis_projects/8363170_Yang/comparisons.csv \
  -s "deseq2 edger noiseq"

# p. aeruginosa

sbatch /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/8363170_Yang/paeruginosa-differentials-bacterial-short30 \
  -g /data/gencore/analysis_projects/8363170_Yang/paeruginosa-quants-bacterial-short30/gene_count_matrix.csv \
  -c /data/gencore/analysis_projects/8363170_Yang/comparisons.csv \
  -s "deseq2 edger noiseq"
```

The next step is to merge the results into a summary output of the three tools, and create some Venn diagrams to visualize the overlap in DEG calls.

```
source activate /data/biocore/programs/mamba-envs/biocore-rna
python /data/gencore/shared_scripts/github-repos/project-scripts/Referenced_Scripts/RNA-DEG_Modules_2025/Module_B/merge_de_both.py /data/gencore/analysis_projects/8363170_Yang/comparisons.csv 0
```

## Module C: Functional Enrichment

Because these are non-standard genomes (or at least they aren't easily accessible in the programs I use for GO comparisons), I use the GFF files to create a GO to Genes correlation table for the genome.

```
# e. coli
grep 'GO' CFT073.genomic.gff | grep -o '\Parent[^;]*' | awk -F '-' '{ print $2 }' > gene-ids.txt
grep -o '\Ontology_term[^;]*' CFT073.genomic.gff | awk -F '=' '{ print $2 }' > gos.txt

paste gene-ids.txt gos.txt > CFT073.ids-and-gos.txt

awk 'BEGIN{FS=OFS="\t"}
  {
    gos=$2;
    len=length(gos);
    for(i=1; i<len; i=i+11){
      $2="\"" substr(gos, i, 10) "\"";
      print;
    }
  }' CFT073.ids-and-gos.txt > CFT073.ids-and-gos.oneline.txt

sort CFT073.ids-and-gos.oneline.txt | uniq | awk -v OFS='\t' '{ print $2, $1 }' | sed 's/\"//g' > CFT073.goterms.txt
```

So, after I did this I ran into some errors trying to use my custom scripts from before. Looking back at another experiment, I saw a simpler way to create this file that should hopefully avoid those errors! The gene ontology GAF file can be downloaded via the FTP option on the genome page at NCBI.

```
grep 'GO' GCF_014262945.1_ASM1426294v1_gene_ontology.gaf | awk '{print $5,$3}' | tr ' ' '\t' > CFT073-fast-goterms.txt
```


None of the annotation files for P. aeruginosa contained GO terms, so I had the option either to run Interproscan to identify some or just omit functional profiling. Since it's a bacterial genome and Jiseon is a long-time collaborator, I chose to run Interproscan.

```
# p.aeruginosa

#!/bin/bash

##### interproscan functional annotation of genome assemblies #####

#SBATCH -o slurm.%j.interproscan.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.interproscan.err               # STDERR (%j = JobId)
#SBATCH -p htc
#SBATCH -q public
#SBATCH -t 0-4:00                    # estimated time needed - this is quite fast since this is a small genome
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G

module load openjdk-17.0.3_7-gcc-12.1.0 # sol

cd /data/gencore/databases/reference_genomes/pseudomonas_aeruginosa

faa=/data/gencore/databases/reference_genomes/pseudomonas_aeruginosa/GCF_000006765.1_PAO1.protein.faa

echo "scanning sample P. aeruginosa genome in file $faa"
  /data/biocore/programs/interproscan-5.69-101.0/interproscan.sh -i "$faa" -t p \
                                                      -iprlookup -goterms \
                                                      -verbose -b ./GCF_000006765.1_PAO1.ips-out
```

```
awk -F "\t" '{ print $9 }' GCF_000006765.1_PAO1.genomic.gtf | grep -o '\gene_id [^;]*' | awk '{ print $2 }' > gene_ids.txt

awk -F "\t" '{ print $9 }' GCF_000006765.1_PAO1.genomic.gtf | grep -o '\protein_id [^;]*' | awk '{ print $2 }' > protein_ids.txt

paste gene_ids.txt protein_ids.txt |sort | uniq > gene-to-protein.txt

sed -i 's/"//g' gene-to-protein.txt

grep -o '\Ontology_term[^;]*' GCF_000006765.1_PAO1.ips-out.gff3 | awk -F '=' '{ print $2 }' > ips-gos.txts

grep '\Ontology_term[^;]*' GCF_000006765.1_PAO1.ips-out.gff3 | awk '{print $1}' > ips-protein-ids.txt

paste ips-protein-ids.txt ips-gos.txts > protein-ids-to-gos.txt

awk 'BEGIN{FS=OFS="\t"}
  {
    gos=$2;
    len=length(gos);
    for(i=2; i<len; i=i+13){
      $2="\"" substr(gos, i, 10) "\"";
      print;
    }
  }' protein-ids-to-gos.txt > protein-ids-to-gos-oneline.txt

sort protein-ids-to-gos-oneline.txt | uniq > protein-ids-to-gos-nodups.txt

sed -i 's/"//g' protein-ids-to-gos-nodups.txt

join -1 2 -2 1 -o 1.1,2.2 <(sort -k2 gene-to-protein.txt) <(sort -k1 protein-ids-to-gos-nodups.txt) > gene-to-go.txt

join -1 2 -2 1 -o 2.2,1.1 <(sort -k2 gene-to-protein.txt) <(sort -k1 protein-ids-to-gos-nodups.txt) > go-to-gene.txt

sed -i 's/ /\t/g' go-to-gene.txt

cp go-to-gene.txt PAO1.goterms.txt
```


```
grep 'GO' CFT073.genomic.gff | grep -o '\Parent[^;]*' | awk -F '-' '{ print $2 }' > gene-ids.txt
grep -o '\Ontology_term[^;]*' CFT073.genomic.gff | awk -F '=' '{ print $2 }' > gos.txt

paste gos.txt gene-ids.txt > CFT073.ids-and-gos.txt

awk 'BEGIN{FS=OFS="\t"}
  {
    gos=$2;
    len=length(gos);
    for(i=1; i<len; i=i+11){
      $2="\"" substr(gos, i, 10) "\"";
      print;
    }
  }' CFT073.ids-and-gos.txt > CFT073.ids-and-gos.oneline.txt
sort CFT073.ids-and-gos.oneline.txt | uniq > CFT073.goterms.txt
```

Once the GO terms are established, I can proceed with the functional annotation. Those scripts are included in separate R files.
