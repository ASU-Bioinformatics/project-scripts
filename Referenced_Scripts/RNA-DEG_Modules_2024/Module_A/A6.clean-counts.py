# takes gene_count_matrix as input and outputs two new files
# first, a gene_count_matrix with only IDs in the column (cleans the gene_id data)
# second, a list of IDs with corresponding common names

# import necessary python packages
import os, sys
from pathlib import Path
import argparse
import math
import csv

# establish arguments
parser = argparse.ArgumentParser()

parser.add_argument("m", type=str,
                    help="file name of gene count matrix")

args = parser.parse_args()

print("reading in gene count matrix")

mtxF = open(args.m,'r')
mtxData = mtxF.readlines()

counts = []
names = []

for line in mtxData:
    if "," in line:
        newline = line.split(",")
        newline[len(newline)-1]=newline[len(newline)-1].rstrip("\n%")
        totalID = newline[0].split("|")
        geneID = totalID[0].split(".")[0]
        if len(totalID) > 1:
            geneName = totalID[1]
        else:
            geneName = "NA"
        geneLine = newline[1:]
        geneLine.insert(0, geneID)
        print(geneLine)
        counts.append(geneLine)
        names.append(totalID)


with open("gene_count_matrix_idsonly.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(counts)

with open("gene_names.csv", 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerows(names)
