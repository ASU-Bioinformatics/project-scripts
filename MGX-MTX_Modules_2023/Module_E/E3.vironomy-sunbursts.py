#this script takes the level-7 tsv file output by Qiime2/Biom as input
#it outputs a csv for each sample in the analysis containing taxonomy and counts as a csv file
#it also outputs an html file for each sample containing that data as a sunburst plot
# the environment /data/biocore/programs/conda-envs/trapp-env has all necessary packages on Sol

import os, sys
from pathlib import Path
import csv
import pandas as pd
import plotly.graph_objects as go
import plotly.express as px
import argparse
import functools
import textwrap

# establish arguments
parser = argparse.ArgumentParser()

parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter,
                                 usage = "sunburst_csvs.py /path/to/input",
                                 description = textwrap.dedent('''\
                                    Example:
                                    python sunburst_csvs.py /data/gencore/vironomy/table.csv

                                    give an absolute path to the prepared data table
                                    '''))

parser.add_argument("pathToFile", type=str)

args = parser.parse_args()

input = args.pathToFile

charSplit = ";"

taxaF = open(input)
taxaList = taxaF.readlines()
taxaF.close()
taxaArray = []
for line in taxaList:
    taxaArray.append(line.rstrip().split(","))

prefix_dict = {
    0: ["p__"],
    1: ["c__"],
    2: ["o__"],
    3: ["f__"],
    4: ["g__"]
}

for contig in range(1,len(taxaArray[0])):
    print(taxaArray[0][contig])
    count_data = []
    for x in range(1, len(taxaArray)):
        count_data.append([])
        place = len(count_data)-1
        taxon = taxaArray[x][0].rstrip().split(charSplit)
        for t in range(0,5):
            entry = taxon[t]
            for p in prefix_dict[t]:
                entry = entry.removeprefix(p)
                count_data[place].append(entry)
        count_data[place].append(taxaArray[x][contig].rstrip("\n"))
    with open(taxaArray[0][contig]+".csv", 'w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(count_data)
    df = pd.DataFrame(data=count_data, index=None, columns=["phylum", "class", "order", "family", "genus", "abundance"])
    print(df)
    df["abundance"] = pd.to_numeric(df["abundance"])
    df = (df.groupby(["phylum", "class", "order", "family", "genus"], sort=False, as_index=False, dropna=False)
            .agg({'abundance':'sum'})
            .reindex(columns=df.columns))

    score = df['abundance'].sum()/len(df['abundance'])
    figX = go.Figure()
    figX = px.sunburst(df,
                    path=["phylum", "class", "order", "family", "genus"],
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
    figF.write_html(taxaArray[0][contig]+".html",include_plotlyjs='directory')
