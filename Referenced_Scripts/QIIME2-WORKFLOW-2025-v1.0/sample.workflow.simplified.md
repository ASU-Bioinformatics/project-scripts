## ITS Amplicon Analysis of Microbial Populations Using Qiime2

This workflow file provides step-by-step informatics methods. I've attempted to make the code easy to follow, but note that the output structure of some of these scripts, especially the Qiime2 scripts, doesn't match the results folder structure I returned in the SFTP folder; the scripts are designed for computational utility while the return folder is designed for human readability.  

### Cut Adapters and Trim/Filter by Quality

The goal of this script is to remove any adapter sequences from either end of the reads, which should decrease the number of reads that DADA2 flags as chimeric and removes from the downstream analysis. This is especially valuable for paired-end DADA2, I think.

Since this particular run was 2x250bp I used the default cut length of 250bp. The minimum read length is also the default, 50. Other parameter details can be found inside the script itself.

The files output to the cutadapt and unpaired fastq files are intermediates that were not returned or used downstream.

```
sbatch ./scripts/cut-trim-filter_args.sh \
        -i ./raw-fastq \
        -c /tmp/cut-fastq \
        -p ./trimmed-fastq \
        -u /tmp/unpaired-fastq \
        -a ./scripts/adapters.fa
```

### Check QC

After trimming off the adapters, it's good to run fastqc and multiqc. I like to run the QC checks on the original fastq as well, if this hasn't already been done, to have a baseline for comparison. This gives valuable information about the remaining read counts, lengths, and quality, and verifies read depth before continuing the analysis.

```
sbatch ./scripts/fastqc_generation_args_phx.sh \
        -f ./trimmed-fastq

sbatch ./scripts/fastqc_generation_args_phx.sh \
        -f ./raw-fastq
```

### Qiime2 Part 1: Denoising and Metrics

For this step, I imported the paired-end fastq reads into Qiime2 and carried out denoising using paired-end DADA2. If amplicon length had prevented DADA2 from merging the forward and reverse reads into longer amplicons for downstream analysis for a large percentage of reads, I would have used single-end DADA2 on just the forward reads; however, when possible paired-end DADA2 is more accurate.

```
sbatch ./scripts/qiime2_pt1_args.sh \
          -f ./trimmed-fastq \
          -q ./qiime-stats \
          -m ./scripts/Tyson_Terry_Metadata.txt \
          -p "p" -o "demux dada2 stats"
```

Fortunately, with a minimum of 480,058 reads passing DADA2 filters for each sample, there were far more than enough reads to use the paired-end denoising and merging.

### Qiime2 Part 2: Diversity and Taxonomy

#### Subsection 1

This script creates all the alpha and beta diversity metrics, as well as the taxonomic classifications and bar plots.

The value for `-x` is the smallest feature count for any of the samples in the set, obtained from the table QZV file created in the previous step. I specified "Treatment" as the discrete category to analyze and "Above-160F" and "Total" as the continuous categories.

```
sbatch qiime2_pt2_args.sh \
          -q ./qiime-pairedDADA2 \
          -r ./scripts/unite_ver10_99_s_all_04.04.2024-Q2-2024.5.qza \
          -m ./scripts/Tyson_Terry_Metadata.txt \
          -c "Treatment" -n "Above-160F Total" \
          -s "p" -x 480058 -d 100
```

#### Subsection 2

Next, we performed ANCOM-BC analysis using the taxonomic classifications to determine differential abundances between groups in different categories. This is a more manual script so the project-specific variables are included in the code itself, at `qiime2_ancom_args.sh` in the scripts folder. As the only discrete category, only the treatment groups could be compared this way.

ANCOM-BC results were calculated for all seven taxonomic levels, and both visual and tabular data are provided in the return folder.

#### Subsection 3

We also have included interactive sunburst plots to quickly visualize the taxonomic composition of each sample at different hierarchical levels. These were created using the included python script `sunburst_csvs_v2.py`.

```
python ./scripts/sunburst_csvs_v2.py \
  ./taxonomy-analyses/taxonomic-classifications/rel-table6-paired-paired \
  rel-level6-table-paired-paired.tsv
```

As long as the `csvs` and `htmls` folders remain in the same directory and aren't renamed, and html file can be opened in your browser. Changing the names or separating the two folders will interfere with the json script that allows the plots to be viewable and interactive.

### Additional Notes

#### Citations for Published Tools

##### FastQC
Andrews S, Krueger F, Segonds-Bichon A, Biggins L, Krueger C, Wingett S. FastQC: a quality control tool for high throughput sequence data. Babraham Institute. Published online 2010. http://www.bioinformatics.babraham.ac.uk/projects/fastqc

##### MultiQC
Ewels P, Magnusson M, Lundin S, Käller M. MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics. 2016;32(19):3047-3048. doi:10.1093/bioinformatics/btw354

##### QIIME2
Bolyen E, Rideout JR, Dillon MR, et al. Reproducible, interactive, scalable and extensible microbiome data science using QIIME 2. Nat Biotechnol. 2019;37(8):852-857. doi:10.1038/s41587-019-0209-9

##### UNITE Fungal Classifier
Abarenkov, Kessy; Zirk, Allan; Piirmann, Timo; Pöhönen, Raivo; Ivanov, Filipp; Nilsson, R. Henrik; Kõljalg, Urmas (2025): UNITE QIIME release for Fungi. Version 19.02.2025. UNITE Community. doi:10.15156/BIO/3301242

#### Questions?

Please contact the core at asubioinformatics@asu.edu or kristina.buss@asu.edu if you have any questions or concerns about methods or results!
