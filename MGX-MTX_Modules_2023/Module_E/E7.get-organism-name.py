# this script takes a list of genbank ids and retrieves corresponding annotations from Entrez in tabular format
# specifically, it retrieves the organism name, the molecule type, the topology, and the taxonomy
# (molecule type and toplogy are included since this is primarily intended for viral annotations)

# run with list of genbank ids from vironomy's secondary clustering csv file
# can be run in the nanopore-env mamba environment on /data/biocore
# the only argument needed is the full pathname of the file containing the list of ids.

import os, sys
from pathlib import Path
import csv
from Bio import Entrez
from Bio import SeqIO
from urllib.error import HTTPError

Entrez.email = 'kristina.buss@asu.edu'
Entrez.api_key = '88f1bc3fefb8bfba181c831a5def97828309'

#print(sys.argv)
#idFile = open(sys.argv[1], "a")
#idList = idFile.readlines()
#idFile.close()

with open(sys.argv[1], "r") as file:
    idList = file.readlines()

annotArray = []
f = []

for orgID in idList :
    #print(orgID)
    orgID = orgID.rstrip()
    try:
        handle = Entrez.efetch(db="nucleotide", id=orgID, rettype="gb", retmode="text")
    except HTTPError:
        print("HTTPError raised for ", orgID)
    else:
        x = SeqIO.read(handle, 'genbank')
        f.append(x)
        #print(x)
        newLine = [orgID, x.annotations['organism'], x.annotations['molecule_type'], x.annotations['topology'], x.annotations['taxonomy']]
        annotArray.append(newLine)

with open("entrez-output.csv", 'w', newline='') as newFile:
    writer = csv.writer(newFile)
    writer.writerows(annotArray)

with open("entrez-fasta.fa", "w") as output_handle:
    SeqIO.write(f, output_handle, "fasta")

print(annotArray[0:10])
