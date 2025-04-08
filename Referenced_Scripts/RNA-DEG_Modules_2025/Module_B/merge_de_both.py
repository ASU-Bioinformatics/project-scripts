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
                                 usage = "merge_de.py comparisonFile mode",
                                 description = textwrap.dedent('''\
                                    Example:
                                    python merge_de.py comparisons.csv foldChangeCutoff

                                    the comparison file must be in csv format, with a header row containing the names of the comparisons
                                    comparison names should match the file names of the DEG lists

                                    The fold change cutoff should be a number!

                                    '''))

parser.add_argument("comparisonFile", type=str,
                    help="csv file containing the list of comparisons used for DEG analysis")

parser.add_argument("foldChangeCutoff", type=float)

args = parser.parse_args()


compF = open(args.comparisonFile)
compDat = compF.readlines()
compF.close()
compList = compDat.pop(0)
compSplit = compList.split(",")
compSplit.pop(0)

log2fc = args.foldChangeCutoff

# these are default for unannotated files
ePadjField = 5
eLogField = 1
nPadjField = 4
nLogField = 5
dPadjField = 6
dLogField = 2

print(compList)

conditions = ("up", "down")

for condition in conditions:
    for comp in compSplit:
        print(comp)
        comp = comp.rstrip()
        edgeRF = open(comp + ".deg." + condition + ".sig.edgeR.txt")
        edgeRList = edgeRF.readlines()
        edgeRF.close()
        edgeRList.pop(0)
        eList = []
        eDictPadj = {}
        eDictLogFC = {}
        if (edgeRList != ['\n']):
            for line in edgeRList:
                esplit = line.rstrip().split("\t")
                if abs(float(esplit[eLogField])) >= log2fc :
                    eList.append(esplit[0])
                    eDictPadj[esplit[0]] = esplit[ePadjField]
                    eDictLogFC[esplit[0]] = esplit[eLogField]
                else:
                    #print("edgeR log2fc too small: " + str(float(esplit[1])))
                    continue

        noiseqF = open(comp + ".deg." + condition + ".sig.noiseq.txt")
        noiseqList = noiseqF.readlines()
        noiseqF.close()
        noiseqList.pop(0)
        nList = []
        nDictLogFC = {}
        nDictPadj = {}
        if (noiseqList != ['\n']):
            for line in noiseqList:
                nsplit = line.rstrip().split("\t")
                if abs(float(nsplit[nLogField])) >= log2fc :
                    nList.append(nsplit[0])
                    nDictLogFC[nsplit[0]] = nsplit[nLogField]
                    nDictPadj[nsplit[0]] = str(1 - float(nsplit[nPadjField]))
                else:
                    continue

        deseq2F = open(comp + ".deg." + condition + ".sig.deseq.txt")
        deseq2List = deseq2F.readlines()
        deseq2F.close()
        deseq2List.pop(0)
        dList = []
        dDictPadj = {}
        dDictLogFC = {}
        if (deseq2List != ['\n']):
            for line in deseq2List:
                dsplit = line.rstrip().split("\t")
                if abs(float(dsplit[dLogField])) >= log2fc :
                    dList.append(dsplit[0])
                    dDictLogFC[dsplit[0]] = dsplit[dLogField]
                    dDictPadj[dsplit[0]] = dsplit[dPadjField]
                else:
                    #print("deseq2 log2fc too small: " + str(float(dsplit[6])))
                    continue

        plt.figure(figsize=(6,6))
        venn = venn3([set(eList), set(nList), set(dList)], set_labels = ('edgeR', 'NOISeq', 'DESeq2'))
        plt.title('DE Genes: ' + comp + " " + condition + "-regulated", fontsize = 12)

        plt.savefig("venn." + comp + ".deg." + condition + "." + str(log2fc).split('.')[0] + ".venn.png")
        plt.close()

        unique_DEs = set(eList + nList + dList)
        de_list = [["gene_id","edgeR-padj","NOISeq-padj","DESeq2-padj","average-padj","edgeR-log2fc","NOISeq-log2fc","DESeq2-log2fc","average-log2fc"]]
        for de in unique_DEs:
            de_line = [de,"NA","NA","NA","NA","NA","NA","NA","NA"]
            sumPadj = 0.0
            sumLogFC = 0.0
            number = 0
            if eList.count(de) > 0:
                de_line[1] = eDictPadj[de]
                de_line[5] = eDictLogFC[de]
                sumPadj = sumPadj + float(eDictPadj[de])
                sumLogFC = sumLogFC + float(eDictLogFC[de])
                number = number + 1
            if nList.count(de) > 0:
                de_line[2] = nDictPadj[de]
                de_line[6] = nDictLogFC[de]
                sumPadj = sumPadj + float(nDictPadj[de])
                sumLogFC = sumLogFC + float(nDictLogFC[de])
                number = number + 1
            if dList.count(de) > 0:
                de_line[3] = dDictPadj[de]
                de_line[7] = dDictLogFC[de]
                sumPadj = sumPadj + float(dDictPadj[de])
                sumLogFC = sumLogFC + float(dDictLogFC[de])
                number = number + 1
            avgPadj = sumPadj / number
            avgLogFC = sumLogFC / number
            de_line[4] = avgPadj
            de_line[8] = avgLogFC
            de_list.append(de_line)

        with open(("merged_deg." + comp + "." + condition + "." + str(log2fc) + ".stats.csv"), 'w', newline='') as file:
            writer = csv.writer(file)
            writer.writerows(de_list)
