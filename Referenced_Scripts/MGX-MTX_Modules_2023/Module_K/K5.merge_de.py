#takes DE output from three different analyses and creates a Venn Diagram

import os, sys
import csv
from pathlib import Path
from matplotlib_venn import venn3, venn3_circles
from matplotlib import pyplot as plt
import argparse
import functools
import textwrap

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "merge_de.py comparisonFile deseqDir edgeRDir noiseqDir",
                                 description = textwrap.dedent('''\
                                    Example:
                                    python rna_gene_info.py comparisons.csv ./automated-deseq ./automated-edger ./automated-noiseq

                                    the geneInformation file must be in csv format, with a header row containing the names of the comparisons
                                    comparison names should match the file names of the DEG lists
                                    the three directories contain the output from the three DEG tools
                                    '''))

parser.add_argument("comparisonFile", type=str,
                    help="csv file containing the list of comparisons used for DEG analysis")

parser.add_argument("deseqDir", type=str,
                    help="the directory containing DESeq2 output")

parser.add_argument("edgeRDir", type=str,
                    help="the directory containing edgeR output")

parser.add_argument("noiseqDir", type=str,
                    help="the directory containing NOIseq output")

args = parser.parse_args()

compF = open(args.comparisonFile)
compDat = compF.readlines()
compF.close()
compList = compDat.pop(0)
compSplit = compList.split(",")
compSplit.pop(0)

print(compList)

conditions = ("up", "down")

for condition in conditions:
    for comp in compSplit:
        comp = comp.rstrip()
        edgeRF = open(args.edgeRDir + "/" + comp + ".deg." + condition + ".sig.edgeR.txt")
        edgeRList = edgeRF.readlines()
        edgeRF.close()
        edgeRList.pop(0)
        eList = []
        eDict = {}
        if (edgeRList != ['\n']):
            for line in edgeRList:
                esplit = line.split("\t")
                eList.append(esplit[0])
                eDict[esplit[0]] = esplit[1]

        noiseqF = open(args.noiseqDir + "/" + comp + ".deg." + condition + ".sig.noiseq.txt")
        noiseqList = noiseqF.readlines()
        noiseqF.close()
        noiseqList.pop(0)
        nList = []
        nDict = {}
        if (noiseqList != ['\n']):
            for line in noiseqList:
                nsplit = line.split("\t")
                nList.append(nsplit[0])
                nDict[nsplit[0]] = nsplit[3]

        deseq2F = open(args.deseqDir + "/" + comp + ".deg." + condition + ".sig.deseq.txt")
        deseq2List = deseq2F.readlines()
        deseq2F.close()
        deseq2List.pop(0)
        dList = []
        dDict = {}
        if (deseq2List != ['\n']):
            for line in deseq2List:
                dsplit = line.split("\t")
                dList.append(dsplit[0])
                dDict[dsplit[0]] = dsplit[2]

        plt.figure(figsize=(4,4))
        venn = venn3([set(eList), set(nList), set(dList)], set_labels = ('edgeR', 'NOISeq', 'DESeq2'))
        plt.title('DE Genes: ' + comp + " " + condition + "-regulated")

        plt.savefig(comp + ".deg." + condition + ".png")
        plt.close()

        unique_DEs = set(eList + nList + dList)
        de_list = [["gene_id","edgeR","NOISeq","DESeq2"]]
        for de in unique_DEs:
            de_line = [de,"NA","NA","NA"]
            if eList.count(de) > 0:
                de_line[1] = eDict[de]
            if nList.count(de) > 0:
                de_line[2] = nDict[de]
            if dList.count(de) > 0:
                de_line[3] = dDict[de]
            de_list.append(de_line)

        with open(("merged_deg." + comp + "." + condition + ".csv"), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(de_list)
