# Standard Variant Calling Data Return

The goal of this analysis is to identify SNPs and structural variants within one or more genomes as compared to a reference genome. This can either be a known reference from NCBI or an in-house reference genome; if de novo assembly of an in-house reference is needed, that would be part of a separate analysis request. Ploidy is assumed to be diploid unless otherwise noted in the bioinformatics request.

After these variants are detected, the impact of each variant is predicted and gene-specific annotations are generated for each one. If desired, genotyping and lineage analysis can also be included.

## Alignment Output

[BWA]("https://doi.org/10.48550/arXiv.1303.3997", "Title"), an alignment tool based on the Burrows-Wheeler algorithm, will be used to align your short-read sequences to the reference genome. SAM (human readable) alignment files, sorted BAM (binary compressed) alignment files and BAI indexes for the BAM files are included in the returned data. In addition to these SAM/BAM output files, a flagstats quality metric file is created for each sample. BAM files are useful for visualization,

### How to Read a SAM File

Most of the time, you won't need to open the SAM files - along with the BAM files, they primarily serve as input for downstream analysis or visualization tools. However, it can still be useful to understand their file structure and the meaning of the internal scores and flags.

SAM files are tab-delimited text files with 11 columns required for each read in the input sequencing file for a given sample.

1. QNAME
  * This column contains the name of the input sequencing read, typically the sequence header name assigned automatically by the Illumina sequencing instrument. If two query names are identical, the aligner considers them to be two segments of the same template, ensuring that paired forward and reverse reads are aligned together.
2. FLAG
  * This column is probably the most complex in the alignment file, but also one of the most important for understanding the quality of the alignment. The sum of different binary (bitwise) scores signifying different pieces of information about the alignment, as shown in the following table.

  | Bit Score | Description |
  | --------- | ----------- |
  | 1 | template contains multiple segments (here, both forward and reverse reads are present) |
  | 2 | all segments are properly aligned |
  | 4 | the segment is unmapped (here, the forward read isn't mapped to the reference)|
  | 8 | the next segment in the template is unmapped (here, the reverse read isn't mapped)|
  | 16 | the sequence is being reverse complemented (here, used for the reverse read only) |
  | 32 | the next sequence in the template is being reverse complemented (here, used for the forward read only) |
  | 64 | the first segment in the template (can be either the forward or reverse read) |
  | 128 | the last segment in the template (can be either the forward or reverse read) |
  | 256 | this alignment is a secondary alignment for a sequence that can be mapped to the reference more than once |
  | 512 | this alignment doesn't pass quality filters |
  | 1024 | this sequence is a duplicate from PCR or optical amplification |
  | 2048 | this is a supplementary alignment (the second alignment from a chimeric read) |

  * Let's use a bit score of 99 as an example. The only way to obtain this score with the binary bit score system is as the sum of 64, 32, 2, and 1. Flag 64 tells us this query is the first in a series of segments within one template, flag 32 tells us that the next sequence in the template is being reverse complemented, flag 2 tells us that both reads are properly aligned to the reference, and flag 1 confirms that both reads in the pair are present.
  * Another example is the bit score 147. This score represents the flags 128, 16, 2, and 1. Flag 128 tells us this is the last segment in the query, flag 16 tells us it has been reverse complemented, and flags 2 and 1 confirm that all the segments in the query and present and properly aligned. In the SAM file I used to find these example scores, the score of 99 was from the forward read in a read pair and the score of 147 was from its paired reverse read.
3. RNAME
  * This is the reference sequence name for the alignment. Typically it will be the name of the chromosome to which the query template is aligned.
4. POS
  * This is the start position of the alignment on the reference genome, counting from a value of 1 for the first base in the reference. A POS of 0 indicates an unmapped read.
5. MAPQ
  * This is the mapping quality score. It's calculated as -10 * log<sub>10</sub>(P<sub>w</sub>), where P<sub>w</sub> is the probability that the alignment is incorrect. Given this score, it's possible to work backwards to determine the probability of a incorrect alignment. A score of 60 corresponds to a probability of incorrect alignment of 0.000001, in turn telling us that the probability the alignment is correct is 0.999999 (where a perfect probability score is 1).
6. CIGAR
  * This is a summary of the alignment, using the letters M (alignment match/mismatch), I (insertion to the reference), D (deletion from the reference), S (soft clipping - unaligned sequence at the end of a read, such as adapter sequence), and H (hard clipping - similar to soft clipping except clipped bases aren't included in the sequence fasta in column 10). For example, a CIGAR string of 69S82M means that the first 69 bases are soft-clipped (most likely adapter read-through) and the next 82 are matched (SNPs or exact matches) without insertions or deletions.
7. RNEXT
  * The reference sequence name of the primary alignment of the next read in the query. For a properly aligned pair this should be the same (denoted by '='); for a chimeric pair it could be the name of a different chromosome.
8. PNEXT
  * The position of the primary alignment of the next read in the query. This will be the same as the POS field for the alignment of the next query segment, and in concert with the TLEN value in column 9 can be used to find the length of the gap between template segments.
9. TLEN
  * The observed template length (including all segments and the distance between them); the value is positive for the leftmost section of the template (aligned closer to the beginning of the reference) and negative for the rightmost section of the template (aligned closer to the end of the reference).
10. SEQ
  * The fasta sequence of the query segment used for the alignment.
11. QUAL
  * The ASCII conversion of the phred-33 quality score. The characters A-K represent the highest quality scores (from 99.9% confidence in the base call to 99.99% confidence).

Optional Tags following Column 11:
* NM: the number of mismatches and gaps in the alignment.
* MD: this field describes SNP/indel calling without looking at the entire reference sequence. So, an MD value of 50G49 would signify 50 bases in the query that match the reference, then a G base on the reference that differs from the aligned base in the query, and then 49 bases in the query that match the reference. This can get complex for more and longer insertions and deletions.
* MC: the cigar string for the paired or next segment in the template
* AS: the alignment score generated by the alignment tool
* XS: the score of the best secondary alignment associated with the query template - if this is close to the AS value, multiple valid alignments are being identified.

## How to Read the Flagstat Metrics Files

All lines are expressed as the reads passing the QC filter plus the reads failing QC

* Total number of forward and reverse reads for the sample
* Total number of primary alignments (neither secondary or supplementary, but not necessarily successfully mapped to the reference)
* Total number of secondary alignments
* Total number of supplementary alignments
* Total number of duplicate reads
* Total number of duplicate reads that are also primary alignments
* Total number of reads successfully mapped to the reference (percentage calculated out of the total number of reads)
* Total number of primary alignments successfully mapped to the reference (percentage calculated out of the total number of primary alignments)
* Number of primary alignments with pairs in the sequencing file
* Number of forward primary alignment reads
* Number of reverse primary alignment reads
* Number of properly paired alignments (percentage calculated out of the total number of paired reads)
* Number of reads where both the read and its mate are successfully mapped to the reference
* Number of reads where only one read in the pair is not mapped to the reference (percentage calculated out of the number of reads with pairs in the sequencing file)
* Number of reads whose mate/pair is mapped to a different chromosome
* Number of reads whose mate/pair is mapped to a different chromosome with a mapping quality greater than or equal to 5

I generally rely on the percent of reads mapped to the reference as a useful fast way to determine the success and quality of the alignment.

## Small Variant Calling

Small variant calling - variants less than 50bp long, typically - are detected using GATK's HaplotypeCaller program. Unless otherwise informed by the researcher, the core presumes a diploid genome for the variant calling steps.

### Individual Sample Variants Only
For this analysis, both unfiltered and quality-filtered VCF files are returned. Default filtering parameters for the VCF file are derived from this Evolution and Genomics [tutorial]("http://evomics.org/learning/population-and-speciation-genomics/2020-population-and-speciation-genomics/first-steps-in-genomic-data-analysis/", "Title") and *Identification of Substrain-Specific Mutations by Massively Parallel Whole-Genome Resequencing of Synechocystis sp. PCC6803* by [Kanesaki et al., 2011]("https://doi-org.ezproxy1.lib.asu.edu/10.1093/dnares/dsr042", "Title"). This type of variant calling is appropriate for single sample comparison to a genome, and can't be used for downstream genotyping.

### Genotyped Variant Calls
For this analysis, unfiltered individual VCF files are returned, along with a merged VCF file that genotypes each variant across the sample set (known as a GVCF file). Merging and genotyping are carried out with CombineGVCFs and GenotypeGVCFs from GATK.

### How to Read a VCF File
VCF files can be large and unwieldy, but it is possible to open them in a text editor such as Microsoft Excel to manually inspect variant locations, characteristics, and quality. Every VCF file has a header section detailing the commands given to the generating software and containing descriptions of the metrics provided for each variant. After the header is a table with 9 columns of general information followed by one column for each sample containing the sample-specific metrics.

1. CHROM, the first column, is just the name of the reference chromosome.
2. POS provides the position on the reference where the variant is (or begins, for an insertion).
3. ID could be a name for the SNP variant, but is typically left unspecified with the '.' character.
4. REF provides the sequence for the reference allele.
5. ALT provides the sequence(s) for the variant allele(s).
6. QUAL is the quality score for the variant call in an individual sample VCF, or of the genotype call for a genotyped GVCF.
7. FILTER specifies whether a variant passed or failed the quality filter; if no filter was used, the '.' is used here.
8. INFO contains a lot of information about the variant, which differs between genotyped and individual samples. The descriptions for each abbreviated code are included in the VCF header; they include metrics such as AF for the frequency of each alternate allele and SOR to estimate strand bias in the call.
9. FORMAT explains the meaning and order of the numeric data included for each sample in the following columns. A common FORMAT column in a genotyped VCF file is GT:AD:DP:GQ:PL, where:
  * GT is the called genotype for the sample at that variant
  * AD is the allele read depths for the reference and alternate alleles, in the order listed in the columns REF and ALT.
  * DP is the approximate total read depth at the specified location.
  * GQ is the genotype quality score.
  * PL is the normalized phred-scaled likelihood for the called genotypes, representing the probability that the genotype call is incorrect.
  * As an example of interpreting this data for a specific sample, consider the sample value 0:237,0:237:99:0,1800. GT is 0 - this sample is assigned to the reference genotype. AD is 237 reads for the reference allele and 0 for the alternate allele. DP tells us that the total read depth is 237. GQ gives the genotype quality score of 99. PL tells us the likelihood of this allele not matching the reference is 0, while the likelihood of this allele not matching the alternate is very high at 1800.

## Mid-Range Variant Calling

These variants cannot be integrated into the GATK genotyping pipeline at this time, so are currently only provided for projects where genotyping is not requested (we're working on this!).

[GRIDSS2]("https://github.com/PapenfussLab/gridss", "Title") focuses on slightly longer variants than HaplotypeCaller, in my experience mostly identifying 50-300bp long breakpoints, insertions, and deletions. In addition, it incorporates a default quality filter for its variant calls, so the script outputs both the complete VCF (including poor quality calls) and the filtered VCF (including only calls that pass the GRIDSS2 filter). Detail of the output VCF file format is provided in the previous section.

## Long Structural Variant Identification

These variants cannot be integrated into the GATK genotyping pipeline at this time, so are currently only provided for projects where genotyping is not requested (we're working on this!). If you anticipate or know in advance that a long variant is likely important for your samples, we can search for these variants individually, however.

This script uses the program [DELLY]("https://github.com/dellytools/delly?tab=readme-ov-file", "Title") to identify long structural variants such as breakpoints and major insertions/deletions. Like GRIDSS2, DELLY incorporates its own quality filter so both unfiltered and filtered VCF files are output. Please refer to the short variant calling section for details on the VCF file format.

## snpEff Annotation

We use [snpEff]("https://pcingola.github.io/SnpEff/", "Title") to predict the impact each variant may have on phenotype, based on the type of variant and genomic region in which it occurs (for example, a stop codon in the middle of a coding region will probably have a higher impact than a point mutation in an intron). snpEff output is also returned in VCF format, and the additional functional annotation is included in the INFO column. In instances with high numbers of variants, it can be useful to filter out only high impact variants; for any dataset, it can also be helpful to convert the VCF file into a simple text table for ease of manual examination. If a specific table format is desired and you're not comfortable doing a conversion yourself, let us know the table specs and we will return a table for each VCF as well.
