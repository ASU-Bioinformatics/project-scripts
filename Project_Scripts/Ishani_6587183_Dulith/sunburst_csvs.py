#this script takes the level-7 tsv file output by Qiime2/Biom as input
#it outputs a csv for each sample in the analysis containing taxonomy and counts as a csv file
#it also outputs an html file for each sample containing that data as a sunburst plot

import os, sys
from pathlib import Path
import csv
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px


clsf = sys.argv[1]
wd = sys.argv[2]+"/"+clsf+"-rel-table7/"

charSplit = ";"

sF = open(wd+clsf+"-rel-level7-table.tsv")
sList = sF.readlines()
sF.close()
sArray = []
for line in sList:
    sArray.append(line.rstrip().split("\t"))

prefix_dict = {
    0: ["k__","d__"],
    1: ["p__"],
    2: ["c__"],
    3: ["o__"],
    4: ["f__"],
    5: ["g__"],
    6: ["s__"]
}

os.makedirs(wd+"csvs/")
os.makedirs(wd+"htmls/")

for s in range(1,len(sArray[1])):
    print(sArray[1][s])
    count_data = []
    for x in range(3, len(sArray)):
        if sArray[x][s] != "0.0":
            count_data.append([])
            place = len(count_data)-1
            taxon = sArray[x][0].rstrip().split(charSplit)
            prev=""
            for t in range(0,7):
                entry = taxon[t]
                for p in prefix_dict[t]:
                    entry = entry.removeprefix(p)
                if (entry == "") or (entry == "__"):
                    if (prev == 'other') or (prev == None):
                        count_data[place].append(None)
                        prev=None
                    else:
                        count_data[place].append('other')
                        prev='other'
                else:
                    count_data[place].append(entry)
                    prev = entry
            count_data[place].append(sArray[x][s].rstrip("\n"))
    with open(wd+"csvs/"+sArray[1][s]+".csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(count_data)
    df = pd.DataFrame(data=count_data, index=None, columns=["kingdom", "phylum", "class", "order", "family", "genus", "species", "abundance"])
    df["abundance"] = pd.to_numeric(df["abundance"])
    df = (df.groupby(["kingdom", "phylum", "class", "order", "family", "genus", "species"], sort=False, as_index=False, dropna=False)
            .agg({'abundance':'sum'})
            .reindex(columns=df.columns))

    score = df['abundance'].sum()/len(df['abundance'])
    figX = go.Figure()
    figX = px.sunburst(df,
                    path=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
                    values='abundance',
                    branchvalues='total')
    figF = go.Figure(go.Sunburst(
                    labels=figX['data'][0]['labels'].tolist(),
                    parents=figX['data'][0]['parents'].tolist(),
                    values=figX['data'][0]['values'].tolist(),
                    ids=figX['data'][0]['ids'].tolist(),
                    insidetextorientation='radial',
                    branchvalues='total',
                    marker=dict(
                        colorscale='blues',
                        cmax=0.1,
                        cmin=0.0),
                    hoverinfo='label+current path+value+percent parent'
    ))
    figF.write_html(wd+"htmls/"+sArray[1][s]+".html",include_plotlyjs='directory')
