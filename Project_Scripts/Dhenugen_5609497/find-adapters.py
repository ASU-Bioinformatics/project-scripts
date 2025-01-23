from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
import os, sys
from pathlib import Path

adaptor = "CCGTGACCTCCTACTCACGGCTAACC"
keep = []

for seq_record in SeqIO.parse(sys.argv[1], "fasta-2line"):
    seq = seq_record.seq
    pend = len(seq_record.seq)
    pstart = pend - 26
    minCount = 0
    perm = seq[pstart:pend]
    alignment = pairwise2.align.globalxx(perm, adaptor, score_only=True)
    if alignment >= 19:
        keep.append(seq_record)

SeqIO.write(keep, "dhenugen-filtered.fasta", "fasta-2line")
