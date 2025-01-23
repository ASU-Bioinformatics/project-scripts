# this takes a directory as input and uses all the final.out alignment logs in that directory
# to generate graphical representations of the alignment statistics for the dataset.

# import necessary python packages
import os, sys
from pathlib import Path
import argparse
import pandas as pd
import matplotlib
import math
import re
import csv
import numpy as np

# establish arguments
parser = argparse.ArgumentParser()

parser.add_argument("d", type=str,
                    help="directory where final.out star alignment logs are found")

args = parser.parse_args()

print("reading in data files")

starlist = []
cols=["SampleName", "StartJob", "StartMap", "Finish", "MapSpeed", "InputReads","AvgInputLength",
      "UniquelyMappedReadsNumber", "UniquelyMappedReadsPercent", "AvgMappedLength",
      "TotalSplices", "AnnotatedSplices", "GT/AG_Splices", "GC/AG_Splices", "AT/AC_Splices", "NoncanonicalSplices",
      "MismatchPerBasePercent", "DeletionPerBasePercent", "AvgDeletionLength", "InsertionPerBasePercent", "AvgInsertionLength",
      "MultimappedReadsNumber", "MultimappedReadsPercent", "OverlyMultimappedReadsNumber", "OverlyMultimappedReadsPercent",
      "UnmappedMismatchNumber", "UnmappedMismatchPercent", "UnmappedShortNumber", "UnmappedShortPercent",
      "UnmappedOtherNumber", "UnmappedOtherPercent", "ChimericNumber", "ChimericPercent"]
starRows=[]
for root, dirs, files in os.walk(args.d):
    for file in files:
        filename, extension = os.path.splitext(file)
        if extension == '.out':
            if 'final' in filename:
                starF = open((args.d + "/" + file),'r')
                stardata = starF.readlines()
                newstardata = [re.split("_STAR", filename)[0]]
                starRows.append(re.split("_STAR", filename)[0])
                for line in stardata:
                    if "\t" in line:
                        newline = line.split("\t")[1].rstrip("\n%")
                        newstardata.append(newline)
                starlist.append(newstardata)

starArray = np.asarray(starlist)
with open('STAR_Stats.csv', 'w', newline='') as f:
    writer = csv.writer(f)
    writer.writerow(cols)
    writer.writerows(starArray)

print("generating plots")
#starRows = ["A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R"]
starDF = pd.DataFrame(starArray, columns = cols, index = starRows)

for i in ["MapSpeed", "InputReads","AvgInputLength",
      "UniquelyMappedReadsNumber", "UniquelyMappedReadsPercent", "AvgMappedLength",
      "TotalSplices", "AnnotatedSplices", "GT/AG_Splices", "GC/AG_Splices", "AT/AC_Splices", "NoncanonicalSplices",
      "MismatchPerBasePercent", "DeletionPerBasePercent", "AvgDeletionLength", "InsertionPerBasePercent", "AvgInsertionLength",
      "MultimappedReadsNumber", "MultimappedReadsPercent", "OverlyMultimappedReadsNumber", "OverlyMultimappedReadsPercent",
      "UnmappedMismatchNumber", "UnmappedMismatchPercent", "UnmappedShortNumber", "UnmappedShortPercent",
      "UnmappedOtherNumber", "UnmappedOtherPercent", "ChimericNumber", "ChimericPercent"]:
      starDF[i] = starDF[i].astype(float)

starDF.sort_index(axis=0, inplace=True)
nums = starDF[["SampleName", "UniquelyMappedReadsNumber", "MultimappedReadsNumber",
              "OverlyMultimappedReadsNumber", "UnmappedMismatchNumber", "UnmappedShortNumber", "UnmappedOtherNumber"]]

starPlot1 = nums.plot.barh(stacked = True, mark_right = True, figsize=(8,len(starRows)/2))

starPlot1.legend(bbox_to_anchor=(0,1,1,0), loc="lower left", mode="expand", ncol=2)
starPlot1.get_figure().tight_layout()
starPlot1.get_figure().savefig("STAR_QC_Reads.png", bbox_inches = 'tight')

errors = starDF[["MismatchPerBasePercent", "DeletionPerBasePercent", "InsertionPerBasePercent"]]

starPlot2 = errors.plot.bar(stacked = False, mark_right = True, figsize=(len(starRows)/2,4))

starPlot2.legend(bbox_to_anchor=(0,1,1,0), loc="lower left", mode="expand", ncol=2)
starPlot2.get_figure().tight_layout()
starPlot2.get_figure().savefig("STAR_QC_Errors.png", bbox_inches = 'tight')
