sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/A1.alignment.sh \
  -a /data/gencore/analysis_projects/7055035_Misra/bwa-align \
  -r /data/gencore/databases/reference_genomes/ecoli \
  -g HG738867.1.genomic.fasta \
  -f /data/gencore/analysis_projects/7055035_Misra/fastq

sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/B4.gatk-gvcf.sh \
        -i /data/gencore/analysis_projects/7055035_Misra/bwa-align \
        -r /data/gencore/databases/reference_genomes/ecoli/HG738867.1.genomic.fasta \
        -o /data/gencore/analysis_projects/7055035_Misra/variant-calls-haploid \
        -n 'misra.all'

#snp database to use: Escherichia_coli_str_k_12_substr_mc4100

sbatch /data/gencore/shared_scripts/variant_calling/updated-pathway-2023/C1.snpeff-annotation.sh \
  -i /data/gencore/analysis_projects/7055035_Misra/variant-calls-haploid \
  -r Escherichia_coli_str_k_12_substr_mc4100 \
  -o /data/gencore/analysis_projects/7055035_Misra/variant-calls-haploid \
  -s /data/biocore/programs/snpEff
