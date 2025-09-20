# module D is for non-standard functional analysis, when I don't have a genome available in clusterProfiler
# or just can't find one...

# panther input: list of gene ids from the DEG lists, filtered however works best
# sometimes genes identified by two of the three tools, for example

# I haven't looked into Panther's API yet, so for now this is all done manually
# for each gene set, run the three GO analyses and the Panther Pathways analysis
# export the resulting tables
# make sure you note which test and multiple testing correction you're using!

# each of the output tables has a 12-line header which has to be removed
# in order for the R scripts to handle the data
# this is easy with the following code (or an appropriate version for your directory tree)

for i in $(find ./ -type f -name "*.txt")
do
 echo $i
 tail -n +12 $i > "${i%.txt}"-forR.txt
done

# once all the files are formatted correctly, move on to pooled-dotplots.R
