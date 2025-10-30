## Terry ITS Cut/Trim Adapters

The goal of this script is to remove any adapter sequences from either end of the reads, which should decrease the number of reads that DADA2 flags as chimeric and removes from the downstream analysis. This is especially valuable for paired-end DADA2, I think.

While the average ITS region is 550bp long (and I would recommend a cut lenth of 290 on a 300bp run), this particular run was a 2x250 so I'm sticking with the default cut length of 250bp. The minimum read length is also the default, 50.

This script is designed to be run on Phx.

```
sbatch /data/gencore/shared_scripts/qiime_scripts/cut-trim-filter_args.sh \
        -i /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/original \
        -c /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/cutadapt_args \
        -p /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args \
        -u /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/unpaired_args \
        -a /data/gencore/databases/trimmomatic/all.fa
```

After trimming off the adapters, it's good to run fastqc and multiqc again to check on the remaining read counts, lengths, and quality. The only crucial sample set to run this on before continuing the analysis is the paired cut and filtered reads, but I like to run it on the cutadapt and original read sets as well.

```
sbatch /data/gencore/shared_scripts/fastqc_generation_args_phx.sh \
        -f /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args

sbatch /data/gencore/shared_scripts/fastqc_generation_args_phx.sh \
        -f /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/cutadapt_args

sbatch /data/gencore/shared_scripts/fastqc_generation_args_phx.sh \
        -f /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/original
```

## Terry ITS Part 1

I'm going to try running this script with single-end DADA2 and paired DADA2 to compare the output.

Since I am running all the samples in a single directory, I'm not bothering to make a fastq manifest ahead of time.

The cut length (truncation length) removes that number of bases from the end of the read and thus eliminates any reads shorter than that value. The default is 0bp (no truncation) which I think is valid here given the trimming and filtering we've already done to remove poor quality at the end of the reads.

```
sbatch qiime2_pt1_args.sh \
          -f /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args/fastq \
          -q /data/gencore/analysis_projects/8468420_Terry_ITS/qiime-singleDADA2 \
          -m /data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt \
          -p "ps" -o "dada2 stats"
```

```
sbatch qiime2_pt1_args.sh \
          -f /data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args/fastq \
          -q /data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2 \
          -m /data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt \
          -p "p" -o "dada2 stats"
```

## Terry ITS Part 2

The value for -x is the smallest feature count for any of the samples in the set.

```
sbatch qiime2_pt2_args.sh \
          -q /data/gencore/analysis_projects/8468420_Terry_ITS/qiime-singleDADA2 \
          -r /data/biocore/qiime2_classifiers/qiime-2024.5/unite_ver10_99_s_all_04.04.2024-Q2-2024.5.qza \
          -m /data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt \
          -c "Treatment" -n "Above-160F Total" \
          -s "ps" -x 675346 -d 100

sbatch qiime2_pt2_args.sh \
          -q /data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2 \
          -r /data/biocore/qiime2_classifiers/qiime-2024.5/unite_ver10_99_s_all_04.04.2024-Q2-2024.5.qza \
          -m /data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt \
          -c "Treatment" -n "Above-160F Total" \
          -s "p" -x 480058 -d 100
```

## ANCOM BC

```
#!/bin/bash

##### qiime2 pt3: ancom #####

#SBATCH -p general
####SBATCH -q grp_kawoodbu
#SBATCH -q public
#SBATCH -o slurm.%j.out                   # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.err                   # STDERR (%j = JobId)
#SBATCH -t 1-0:00                         # estimated time needed (dada2 can take a while)
#SBATCH --mem=64G

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/qiime2-amplicon-2025.7

fastqDir="/data/gencore/analysis_projects/8468420_Terry_ITS/fastq/paired_args/fastq"
qiimeDir="/data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2"
ancomDir="/data/gencore/analysis_projects/8468420_Terry_ITS/qiime-pairedDADA2-ancombc"

metadata="/data/gencore/analysis_projects/8468420_Terry_ITS/Tyson_Terry_Metadata.txt"
mkdir ancomDir

# sometimes these modules don't work if the --help call isn't run first
# it must be a bug
qiime dada2 --help
qiime composition --help

##### ANCOM with Taxa #####

### Classify, Visualize, and Collapse Table with 2024.10 at Class and Genus level ###
# can add different hierarchical levels!

## Formula and reference level need to be manually specified
# I change the file name to reflect these values
# especially when ancombc needs to be run for multiple categories

for i in 1 2 3 4 5 6 7;
do
  qiime composition ancombc \
    --i-table "$qiimeDir"/level"$i"_table-paired-paired.qza \
    --m-metadata-file "$metadata" \
    --p-formula Treatment \
    --p-reference-levels Treatment::control \
    --o-differentials "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza

  qiime composition da-barplot \
    --i-data "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza \
    --p-effect-size-label 'lfc' \
    --p-feature-id-label 'id' \
    --p-error-label 'se' \
    --p-significance-label 'q_val' \
    --p-significance-threshold 0.05 \
    --p-effect-size-threshold 0.0 \
    --p-level-delimiter ';' \
    --o-visualization "$ancomDir"/dabarplot-ancombc-treatment.control-level"$i"_table-paired-paired_qval.05.qzv

  qiime composition tabulate \
    --i-data "$ancomDir"/ancombc-treatment.control-level"$i"_table-paired-paired.qza \
    --o-visualization "$ancomDir"/table-ancombc-treatment.control-level"$i"_table-paired-paired.qzv

  qiime tools export \
    --input-path "$ancomDir"/ancombc-replicate.snsp-level"$i"-2024.10-sol.qza \
    --output-path "$OUTancomDirPUT_DIR"/ancombc-snsp-level"$i"
done
```
