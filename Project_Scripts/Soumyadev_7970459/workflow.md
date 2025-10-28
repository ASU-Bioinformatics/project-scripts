# Scripts for Project 7970459

This project compared the transcriptomes of three species of bacteria (Arthobacter, Bacillus, and Massilia) grown in nitrogen negative and nitrogen positive environments.

Reference assemblies with Prokka annotations were provided by the researcher.

## Alignment

To compensate for RNA degradation, a lot of different parameters were tested with the samples to obtain a usable number of reads while limiting spurious multimappers caused by allowing shorter aligned segments. Increasing the minimum aligned length didn't decrease the mismatch rate, which is higher than desired. The STAR parameters that proved optimal across the species and conditions are as follows:

```
runThreadN 12
readFilesCommand  gunzip -c
outSAMtype  BAM SortedByCoordinate
outSAMstrandField intronMotif
outFilterIntronMotifs RemoveNoncanonical
outFilterMatchNminOverLread 0
outFilterScoreMinOverLread  0
outFilterMatchNmin  20
peOverlapNbasesMin  0
alignIntronMax 1
outFilterMultimapNmax 2
```

The specific code used to call script A for each species (included in `Referenced_Scripts/RNA-DEG_Modules_2025`):

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-fastq \
  -i /scratch/kawoodbu/7970459-arthrobacter-veryshort \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-bacterial-veryshort-2max.txt \
  -r /data/gencore/databases/reference_genomes/arthrobacter/ \
  -a /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-align-bacterial-veryshort-2max \
  -q /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-stringtie-bacterial-veryshort-2max \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-fastq \
  -i /scratch/kawoodbu/7970459-bacillus-veryshort-2max \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-bacterial-veryshort-2max.txt \
  -r /data/gencore/databases/reference_genomes/bacillus/ \
  -a /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-align-bacterial-veryshort-2max \
  -q /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-stringtie-bacterial-veryshort-2max \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A/A.alignment-wrapper.sh \
  -f /data/gencore/analysis_projects/7970459_Soumyadev/massilia-fastq \
  -i /scratch/kawoodbu/7970459-massilia-veryshort-2max \
  -p /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/supplemental_files/star-params-bacterial-veryshort-2max.txt \
  -r /data/gencore/databases/reference_genomes/massilia/ \
  -a /data/gencore/analysis_projects/7970459_Soumyadev/massilia-align-bacterial-veryshort-2max \
  -q /data/gencore/analysis_projects/7970459_Soumyadev/massilia-stringtie-bacterial-veryshort-2max \
  -s /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_A
```

## DEG Analysis

Three different DEG packages in R were used to identify differentially expressed genes in the three comparisons. For this project, high overlap levels between the three tools were observed, suggesting reliable calls. We specifically used the transcript count matrix for analysis because that was the file in which stringtie captured the gene IDs from the Prokka annotations.

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-degs-bacterial-veryshort-2max \
  -g /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-stringtie-bacterial-veryshort-2max/transcript_count_matrix.csv \
  -c /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-comparisons.csv \
  -s "deseq2 edger noiseq"

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-degs-bacterial-veryshort-2max \
  -g /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-stringtie-bacterial-veryshort-2max/transcript_count_matrix.csv \
  -c /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-comparisons.csv \
  -s "deseq2 edger noiseq"

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/B1.DEG.rscripts.sh \
  -d /data/gencore/analysis_projects/7970459_Soumyadev/massilia-degs-bacterial-veryshort-2max \
  -g /data/gencore/analysis_projects/7970459_Soumyadev/massilia-stringtie-bacterial-veryshort-2max/transcript_count_matrix.csv \
  -c /data/gencore/analysis_projects/7970459_Soumyadev/massilia-comparisons.csv \
  -s "deseq2 edger noiseq"
```

For each species, we merged the information from each package to create a summary file for all genes, all significantly over-expressed genes, and all significantly under-expressed genes. The Python scripts for this step are also in `Referenced_Scripts/RNA-DEG_Modules_2025`. merge_de_both.py creates venn diagrams visualizing the overlap between genes called by each of the three tools, in addition to the merged text files.

```
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_notSig.py /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-comparisons.csv 1.0

python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_notSig.py /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-comparisons.csv 1.0

python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_notSig.py /data/gencore/analysis_projects/7970459_Soumyadev/massilia-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/massilia-comparisons.csv 0.0
python /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_B/merge_de_both.py /data/gencore/analysis_projects/7970459_Soumyadev/massilia-comparisons.csv 1.0
```

## Annotation

To extract basic gene name information from the GFF files for annotating the DEG output files, we used the following code.

```
awk -v FS="\t" '{ print $9 }' Massilia.MET4.gff | awk -v FS='[=;]' -v OFS="\t" '{print $2,$4}' | grep 'peg' > massilia.gene.info.txt
awk -v FS="\t" '{ print $9 }' Arthrobacter_O80.gff | awk -v FS='[=;]' -v OFS="\t" '{print $2,$4}' | grep 'peg' > arthrobacter.gene.info.txt
awk -v FS="\t" '{ print $9 }' Bacillus_O64.gff | awk -v FS='[=;]' -v OFS="\t" '{print $2,$4}' | grep 'peg' > bacillus.gene.info.txt
```

```
sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_C/C2.annotate-DEGs.sh \
  --infoFile /data/gencore/databases/reference_genomes/bacillus/bacillus.gene.info.txt \
  --inputDirectory /data/gencore/analysis_projects/7970459_Soumyadev/bacillus-degs-bacterial-veryshort-2max/deglists

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_C/C2.annotate-DEGs.sh \
  --infoFile /data/gencore/databases/reference_genomes/massilia/massilia.gene.info.txt \
  --inputDirectory /data/gencore/analysis_projects/7970459_Soumyadev/massilia-degs-bacterial-veryshort-2max/deglists

sbatch /data/gencore/shared_scripts/RNAseq/RNA-DEGs-Modular-240314-v1.0.2/Module_C/C2.annotate-DEGs.sh \
  --infoFile /data/gencore/databases/reference_genomes/arthrobacter/arthrobacter.gene.info.txt \
  --inputDirectory /data/gencore/analysis_projects/7970459_Soumyadev/arthrobacter-degs-bacterial-veryshort-2max/deglists
```
