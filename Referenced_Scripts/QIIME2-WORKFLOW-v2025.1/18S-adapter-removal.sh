# Parameters for ITS adapter removal

# the first round of cutadapt removes the adapter sequences at the 5' end of each read caused by the 2-step library creation protocol
# as well as the Illumina adapter sequences at the 3' end of each read caused by fragments shorter than the read length

# the second round of cutadapt removes all reads shorter than 30bp and trims the reads to 150bp
# this trim length is for consistency with your past 18S projects, which were run on a 2x150bp flowcell on the MiSeq

cutadapt \
  --nextseq-trim=12 \
  -g "TCGTCGGCAGCGTCAGATGTGTATAAGAGACAGTATCGCCGTTCGGTACACACCGCCCGTC;anywhere;min_overlap=6" \
  -G "GTCTCGTGGGCTCGGAGATGTGTATAAGAGACAGAGTCAGTCAGCATGATCCTTCTGCAGGTT;anywhere;min_overlap=6" \
  -A "CTGTCTCTTATA;anywhere;min_overlap=10" \
  -a "CTGTCTCTTATA;anywhere;min_overlap=10"

cutadapt \
  -m 30 -l 150
