from Bio.Seq import Seq
from Bio import SeqIO
import os, sys
from pathlib import Path
import csv
import matplotlib.pyplot as plt

permutations = {
    "1": "TTAGGG",
    "2": "TAGGGT",
    "3": "AGGGTT",
    "4": "GGGTTA",
    "5": "GGTTAG",
    "6": "GTTAGG"
}

recordSummary = [["SeqID", "First_Permutation", "Repeat_Num1", "End_Reps1", "End_IS", "IS_Length", "IS_Seq", "First_P_After_IS", "Repeat_Num2", "End_Reps2"]]
repeatHistogram = {"Repeats": "Times_Occurred"}
permutationHistogram = [0,0,0,0,0,0]
repsArray = []
filter = []

for seq_record in SeqIO.parse(sys.argv[1], "fasta-2line"):
    seq = seq_record.seq
    #print(seq_record.id)
    resArray = [["RepNum", "Seq", "P1MM", "P2MM", "P3MM", "P4MM", "P5MM", "P6MM", "MinMM", "Perm"]]
    pend = len(seq_record.seq)
    pstart = pend - 6
    reps = 0
    minCount = 0
    firstPerm = "notfound"
    newRecord = ["", "", "", "", "", "", "", "", "", ""]
    while pstart > 0:
        perm = seq[pstart:pend]
        counts = {"1": 0, "2": 0, "3": 0, "4": 0, "5": 0, "6": 0}
        for x in range(len(perm)):
            if perm[x] != permutations["1"][x]:
                counts["1"]+=1
            if perm[x] != permutations["2"][x]:
                counts["2"]+=1
            if perm[x] != permutations["3"][x]:
                counts["3"]+=1
            if perm[x] != permutations["4"][x]:
                counts["4"]+=1
            if perm[x] != permutations["5"][x]:
                counts["5"]+=1
            if perm[x] != permutations["6"][x]:
                counts["6"]+=1
        minP = min(counts, key=counts.get)
        minNum = counts[minP]
        reps+=1
        newLine = [reps, perm, counts["1"], counts["2"], counts["3"], counts["4"], counts["5"], counts["6"], minNum, permutations[minP]]
        resArray.append(newLine)
        pend = pend - 6
        pstart = pend - 6
        if firstPerm == "notfound":
            if newRecord[0] == "":
                if minNum == 0:
                    firstPerm = int(minP)-1
                    permutationHistogram[firstPerm]+=1
            elif minNum == 0:
                firstPerm = int(minP)-1
                newRecord[4] = pend+6
                reps = 0
        elif minNum >= 3:
            minCount+=1
            if minCount == 3:
                if newRecord[0] == "":
                    repsArray.append(reps-3)
                    if (reps-3) in repeatHistogram:
                        repeatHistogram[(reps-3)]+=1
                    else:
                        repeatHistogram[(reps-3)]=1
                    newRecord=[seq_record.id, permutations[str(firstPerm+1)], (reps-3), (pstart+18), "", "", "", "", "", ""]
                    firstPerm = "notfound"
                    reps = 0
                    minCount = 0
                elif newRecord[4] == "":
                    #print(newRecord)
                    #recordSummary.append(newRecord)
                    break
                else:
                    if reps >= 10:
                        newRecord[7] = permutations[str(firstPerm+1)]
                        newRecord[8] = (reps-3)
                        newRecord[9] = (pstart+18)
                        newRecord[5] = newRecord[3] - newRecord[4]
                        newRecord[6] = seq[newRecord[4]:newRecord[3]]
                        #print(newRecord)
                        #recordSummary.append(newRecord)
                        break
                    else:
                        newRecord[4] = ""
                        #recordSummary.append(newRecord)
                        #print(newRecord)
                        break
        else:
            minCount = 0
    #print(newRecord)
    if newRecord[0] == "":
        newRecord[0] = seq_record.id
    else:
        filter.append(seq_record)
    recordSummary.append(newRecord)
    #with open(("tel-" + seq_record.id + ".csv"), 'w', newline='') as file:
    #    writer = csv.writer(file)
    #    writer.writerows(resArray)

plt.figure()
plt.hist(repsArray)
plt.title("Number of Repeats")
plt.savefig("repeat_histogram.png")

plt.figure()
x = [1,2,3,4,5,6]
y = permutationHistogram
plt.bar(x,y)
plt.title("First Permutation")
plt.savefig("permutation_barplot.png")

with open(("summary.csv"), 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(recordSummary)

SeqIO.write(filter, "telomeres.fasta", "fasta-2line")
