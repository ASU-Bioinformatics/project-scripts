#Somehow this script was deleted. :( :( :(
#I think what it should do is make tsv files containing a list of KEGG terms in each sample
#Using the annotated gene matrix table
#And then feed the tsv files to ko_mapper.py
#To create a multi-sample kegg completeness table
#To go into the kegg-heatmap R script.

# create header for annotated gene count file with sample names
#sed -i.txt '1i contigID\ttranscriptID\tMB-001\tMB-003\tMB-004\tMB-005\tMB-006\tMB-011\tMB-012\tMB-014\tMB-016\tMB-017\tMB-018\tMB-020\tMB-021\tMB-023\tMB-024\tMB-025\tMB-028\tMB-033\tMB-035\tMB-037\tMB-038\tMB-039\tMB-040\tMB-044\tMB-045\tMB-047\tMB-049\tMB-050\tMB-053\tMB-055\tseqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute\tncbiID\tncbiName\tkeggID\tkeggName\ttaxonomy\tgoMF\tgoCC\tgoBP\tinterproscanID\tpfamID\tuniProtEnzyme\tdatabase' annotated.gene.count.matrix.tsv
sed -i.txt '1i contigID\ttranscriptID\t16\t1E\t20\t2D\t3C\t4B\t4\t55\t57\t5A\tP1\tP2\tP5\tP7\tseqname\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattribute\tncbiID\tncbiName\tkeggID\tkeggName\ttaxonomy\tgoMF\tgoCC\tgoBP\tinterproscanID\tpfamID\tuniProtEnzyme\tdatabase' annotated.gene.count.matrix.tsv

declare -A sampleArray

sampleArray+=( ["3"]=MB-001
               ["4"]=MB-003
               ["5"]=MB-004
               ["6"]=MB-005
               ["7"]=MB-006
               ["8"]=MB-011
               ["9"]=MB-012
               ["10"]=MB-014
               ["11"]=MB-016
               ["12"]=MB-017
               ["13"]=MB-018
               ["14"]=MB-020
               ["15"]=MB-021
               ["16"]=MB-023
               ["17"]=MB-024
               ["18"]=MB-025
               ["19"]=MB-028
               ["20"]=MB-033
               ["21"]=MB-035
               ["22"]=MB-037
               ["23"]=MB-038
               ["24"]=MB-039
               ["25"]=MB-040
               ["26"]=MB-044
               ["27"]=MB-045
               ["28"]=MB-047
               ["29"]=MB-049
               ["30"]=MB-050
               ["31"]=MB-053
               ["32"]=MB-055)

sampleArray+=( ["3"]=16_SQP_full
               ["4"]=1E_SQP_full
               ["5"]=20_SQP_full
               ["6"]=2D_SQP_full
               ["7"]=3C_SQP_full
               ["8"]=4B_SQP_full
               ["9"]=4_SQP_full
               ["10"]=55_SQP_full
               ["11"]=57_SQP_full
               ["12"]=5A_SQP_full
               ["13"]=P1_SQP_full
               ["14"]=P2_SQP_full
               ["15"]=P5_SQP_full
               ["16"]=P7_SQP_full )

# the column number is 14 + the sample number
for key in ${!sampleArray[@]}; do
    echo ${key} ${sampleArray[${key}]}
    awk -v col="${key}" -F $'\t' '($col!=0 && $28!="null" && $28!="NA") {print $28}' annotated.gene.count.matrix.tsv > "${sampleArray[${key}]}".keggs.tsv
done

#awk -F $'\t' '($3!=0 && $28!="null" && $28!="NA") {print $28}' annotated.gene.count.matrix.tsv

kolist=$(find ./ -maxdepth 1 -type f -name "*.keggs.tsv" | while read F; do basename $F; done)

python /data/biocore/programs/mamba-envs/microbe-annotator-env/lib/python3.7/site-packages/microbeannotator/pipeline/ko_mapper.py \
  -i $kolist \
  -p ko_map --cluster rows

# find the number of unique Kegg terms identified in at least one sample (the column number is 14 + the number of samples)
awk -F '\t' '{print $28}' annotated.gene.count.matrix.tsv | grep -o -E 'K[[:digit:]]+' | sort | uniq | wc -l
#otak DNA = 8035

# find the number of unique GO MF terms identified in at least one sample (the column number is 17 + the number of samples)
awk -F '\t' '{print $31}' annotated.gene.count.matrix.tsv | grep 'GO' | tr ' ' '\n' | sort | uniq | wc -l
# otak DNA = 4508

# find the number of unique GO CC terms identified in at least one sample (the column number is 18 + the number of samples)
awk -F '\t' '{print $32}' annotated.gene.count.matrix.tsv | grep 'GO' | tr ' ' '\n' | sort | uniq | wc -l
# otak DNA = 1589

# find the number of unique GO BP terms identified in at least one sample (the column number is 19 + the number of samples)
awk '{print $33}' annotated.gene.count.matrix.tsv | grep -o -E 'GO[[:digit:]]+' | tr ' ' '\n' | sort | uniq | wc -l
# otak DNA = 8631
