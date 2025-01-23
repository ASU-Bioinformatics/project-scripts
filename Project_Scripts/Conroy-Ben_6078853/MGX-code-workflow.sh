# E1.vibrant.sh

sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_E/E1.vibrant.sh \
  --contigs /data/gencore/analysis_projects/6078853_Otak_RNA/megahit_alignments/final.contigs.fa \
  --databaseDir /data/gencore/databases/vibrant2/databases \
  --modelDir /data/gencore/databases/vibrant2/files \
  --outDir /data/gencore/analysis_projects/6078853_Otak_RNA/vibrant

# E2.vironomy.sh

#F1 adding missing files
cd /data/gencore/analysis_projects/6078853_Otak_DNA/coassembly-contig-annotations
module load mamba/latest
source activate /data/biocore/programs/mamba-envs/catbat-env
contigs=/data/gencore/sftp/otakuye_conroyben/6078853_Metagenomics/5.metagenome-coassembly/final.contigs.fa
taxaDir=/data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy
databaseDir=/data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database
CAT add_names -i out.CAT.contig2classification.txt \
              -o out.CAT.classification-official.txt \
              -t "$taxaDir" --only_official

CAT summarise -c "$contigs" \
              -i "out.CAT.classification-official.txt" \
              -o out.CAT.summary.txt

CAT contigs -c "$contigs" \
            -d "$databaseDir" \
            -t "$taxaDir" \
            --verbose --index_chunks 1 -n 16 \
            --I_know_what_Im_doing --top 11 \
            -p out.CAT.predicted_proteins.faa \
            -a out.CAT.alignment.diamond

# repeating F module, it looks like it stalled part way through originally
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F1.cat-contig-annotate.sh \
  --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
  --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
  --contigs /data/gencore/sftp/otakuye_conroyben/6078853_Metagenomics/5.metagenome-coassembly/final.contigs.fa \
  --outDir /data/gencore/analysis_projects/6078853_Otak_DNA/contig-redo

# rerun F1 with predicted proteins and alignment from the first time around, according to the log file they are complete
# for this run the -p and -a flags were provided to the CAT function
# I need to find a way to work that into the command line options, if it's important to resume the script.
sbatch /data/gencore/shared_scripts/metagenome_rd/standalone_pathway/Module_F/F1.cat-contig-annotate.sh \
  --databaseDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_CAT_database \
  --taxaDir /data/gencore/databases/CAT_prepare_20210107/2021-01-07_taxonomy \
  --contigs /data/gencore/sftp/otakuye_conroyben/6078853_Metagenomics/5.metagenome-coassembly/final.contigs.fa \
  --outDir /data/gencore/analysis_projects/6078853_Otak_DNA/contig-redo-cat

# determine kegg/go terms from the annotated count matrix
# the column numbers are based on the number of samples, since this is the count matrix. This is for a run with 14 samples.
awk -F "\t" '{ print $28 }' annotated.gene.count.matrix.tsv | sort | uniq > uniq.keggs.txt
wc -l uniq.keggs.txt #9,920, -2 for NA and null makes 9,918

awk -F "\t" '{ print $31 }' annotated.gene.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goMF.txt
wc -l uniq.goMF.txt #4510, -2 for NA and null makes 4,508

awk -F "\t" '{ print $32 }' annotated.gene.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goCC.txt
wc -l uniq.goCC.txt #1591, -2 for NA and null makes 1,589

awk -F "\t" '{ print $33 }' annotated.gene.count.matrix.tsv | tr ' ' '\n' | sort | uniq > uniq.goBP.txt
wc -l uniq.goBP.txt #8633, -2 for NA and null makes 8,631

# to get the taxonomy/other information for the viral alignments (to NCBI refseq database)
# search the list of accession numbers as a text file at: https://www.ncbi.nlm.nih.gov/sites/batchentrez
# then download the results - I think downloading as a GFF3 will be the most successful
# but the download has stalled so far.
