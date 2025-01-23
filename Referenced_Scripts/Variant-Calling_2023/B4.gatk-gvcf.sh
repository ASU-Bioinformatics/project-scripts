#!/bin/bash

##### GATK HaplotypeCaller for genotyping  #####

#SBATCH -p general
#SBATCH -q public
#SBATCH -o slurm.%j.gatkG.out               # STDOUT (%j = JobId)
#SBATCH -e slurm.%j.gatkG.err               # STDERR (%j = JobId)
#SBATCH -t 0-08:00                    # estimated time needed
#SBATCH --mem=64G

#sources for default filtering cutoffs:
#as per recommendations here: https://academic.oup.com/dnaresearch/article/19/1/67/485988?login=true, keep only variants with read depth >10 and >60% of reads representing the alternate sequence
#as per Phred quality score data, keep only variants with GQ >30
#as per GATK recommendations cited here: http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/

module purge

use="DIRECTORY"
outDir="./vcf-output"
ref="./reference.fna"
inputDir="./alignment"
help="FALSE"
params='(FMT/AD[0:1])/(FMT/DP)<0.6 || INFO/DP<10 || FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-10.5 || QD<2.0 || ReadPosRankSum<-8.0'

VALID_ARGS=$(getopt -o hi:r:o:l:p:n: \
                    --long help,inputDir:,ref:,outDir:,list:,parameters:,name: \
                    -- "$@")
if [[ $? -ne 0 ]]; then
  exit 1;
fi

eval set -- "$VALID_ARGS"
while [ : ]; do
  case "$1" in
    -h | --help )
        help="TRUE"
        break
        ;;
    -i | --inputDir)
        echo "Alignment files are found in the directory '$2'"
        inputDir="$2"
        shift 2
        ;;
    -r | --ref)
        echo "The reference genome used is '$2'"
        ref="$2"
        shift 2
        ;;
    -o | --outDir)
        echo "VCF results will be output to the directory '$2'"
        outDir="$2"
        shift 2
        ;;
    -l | --list)
        echo "Alignment files with sample IDs in '$2' will be analyzed"
        use="LIST"
        list="$2"
        shift 2
        ;;
    -n | --name )
        echo "Prefix for joint genotyped output files"
        name="$2"
        shift 2
        ;;
    -p | --parameters)
        echo "The parameter string '$2' will be used to generate the filtered VCF"
        params="$2"
        shift 2
        ;;
    --)
        shift;
        break
        ;;
    *)
        echo "Unexpected option: $1 - please correct."
        ;;
  esac
done

if [ "$help" == "TRUE" ]; then
  cat << EOF
  This script runs GATK HaplotypeCaller with a default sample ploidy of 2 to identify SNPs in WGS samples.
  As input, it requires a sorted BAM file per sample and the reference genome fasta file to which they were aligned
  All samples from the alignment folder will be used unless the -l option provides a quoted, space-separated list of sample ids
  This script is designed for genotype generation/comparison of multiple samples.
  All samples submitted will be jointly genotyped using GATK CombineGCVFs and GenotypeGCVFs

  usage: sbatch B4.gatk-gvcf.sh \
          -i /path/to/input/alignment/files/ \
          -r /path/to/reference/directory/reference.fasta \
          -o /path/to/variant/calling/output \
          -n 'out.prefix.' \
          -p 'INFO/DP<10 || FS>60.0' \
          (-l "sid1 sid2 sid3 ...") (-h)

  options:
    [ -i  |   --inputDir |  directory containing the input alignment files                                                                                ]
    [ -r  |   --ref      |  reference genome in fasta format                                                                                              ]
    [ -o  |   --outDir   |  directory to write output files                                                                                               ]
    [ -n  |   --name     |  name for the combined cohort of samples (will be output prefix)                                                               ]
    [ -p  |   --params   |  parameter string for filtering VCF; it must be single-quoted. The default is:
                            '(FMT/AD[0:1])/(FMT/DP)<0.6 || INFO/DP<10 || FS>60.0 || SOR>3 || MQ<40 || MQRankSum<-10.5 || QD<2.0 || ReadPosRankSum<-8.0]'  ]
    [ -l  |   --list     |   optional; list of sample ids within fastqDir to use for alignment workflow                                                   ]
    [ -h  |   --help     |   prints an informational message and exits script                                                                             ]
EOF
  exit;
fi

module load mamba/latest
source activate /data/biocore/programs/conda-envs/align_breakdance

mkdir -p "$outDir"

if [ "$use" == "DIRECTORY" ];
then
    list=$(find "$inputDir" -maxdepth 1 -type f -name "*.sorted.rdgrp.bam" | \
          while read F; do basename $F; done | cut -d "." -f 1 | sort | uniq)
fi

for i in $list
do
#  echo "Mark Duplicates $i"
#  picard MarkDuplicates I="$inputDir"/"$i".sorted.rdgrp.bam \
#                        O="$inputDir"/"$i".sorted.rdgrp.mrkd.bam \
#                        M=output.metrics CREATE_INDEX=true
#  echo "Index marked bams $i"
#  picard BuildBamIndex INPUT="$inputDir"/"$i".sorted.rdgrp.mrkd.bam
  echo "GATK Haplotype calling $i"
  gatk HaplotypeCaller -R "$ref" \
                      -I "$inputDir"/"$i".sorted.rdgrp.mrkd.bam \
                      -O "$outDir"/"$i".gatk.g.vcf \
                      -A DepthPerAlleleBySample -ERC GVCF \
                      -stand-call-conf 30 --sample-ploidy 1
  echo "$i check coverage"
  bedtools genomecov -ibam "$inputDir"/"$i".sorted.rdgrp.mrkd.bam -max 10 > "$outDir"/"$i".genomecov.txt
done;

# combine gvcf files using CombineGCVFs

cd $outDir
touch variants.list
find "$outDir" -maxdepth 1 -type f -name "*.gatk.g.vcf" > variants.list

gatk CombineGVCFs \
   -R "$ref" \
   --variant variants.list \
   -O "$name".combined.gatk.g.vcf

gatk GenotypeGVCFs \
   -R "$ref" \
   -V "$name".combined.gatk.g.vcf \
   -O "$name".genotyped.gatk.g.vcf

source deactivate
