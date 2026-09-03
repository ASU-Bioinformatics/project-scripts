# *De Novo* Chloroplast Assembly

## Pre-Workflow: Install getOrganelle

This looks like a straightforward installation using Miniconda!

```
module load mamba/latest

mamba create -p /data/biocore/programs/mamba-envs/getOrganelle-env
source activate /data/biocore/programs/mamba-envs/getOrganelle-env

mamba install -c bioconda getorganelle
```

## Step 1: Trim Reads

The initial step is to remove any adapter sequence from the raw reads. The assembler I'm planning to use doesn't recommend running quality control, however. They recommend >5G per end but I'm not sure if that is accounting for any remaining nuclear reads. Hopefully, since this was a chloroplast extraction specifically, there won't be much nuclear contamination.

```
sbatch cut-trim-filter_args.sh \
	--inputDir /data/gencore/analysis_projects/8710018_Sunidhi/fastq \
	--cutadaptDir /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq \
	--adapters /data/gencore/databases/trimmomatic/all.fa

module load fastqc-0.12.1-gcc-11.2.0
cd /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq
fastqc -t 2 *

mkdir -p /data/gencore/analysis_projects/8710018_Sunidhi/cut-qc
cd /data/gencore/analysis_projects/8710018_Sunidhi/cut-qc
mv /data/gencore/analysis_projects/8710018_Sunidhi/cut-fastq/*fastqc* ./

module load mamba/latest
source activate /data/biocore/programs/mamba-envs/multiqc.v1.20/
multiqc *
```

## Step 2: getOrganelle Assemble

A lot of the parameters for this script are automatically detected, while others will need to be optimized. I'm going to run this command within an SBATCH script and check how much time/CPU/memory is consumed.

```
python /data/biocore/programs/mamba-envs/getOrganelle-env/bin/get_organelle_from_reads.py \
	-1 SID_R1.fastq.gz -2 SID_R2.fastq.gz -o SID_chloroplast_output \
	-R 15 -d 21,45,65,85,105 -F emblant_pt
```

I ran multiple settings (included inside the sbatch script) but this setting provided by far the best graph. The getOrganelle selected graph showed the expected large loop / short loop structure with a connecting repeated region, although the repeated region was shorter than typical. There were only a few unresolved loops.

## Step 3: Annotatation

### Annotate with GeSeq

I attempted to annotate with an online tool, GeSeq, which I see referenced in many publications. However, I'm only seeing five of the photosystem I subunits which makes me unsure. The number of predicted genes and the overall length of the sequence are comfortably within the typical range, so I'm wondering if there is just a problem with the gene name annotations by GeSeq. I am trying to run GeSeq again including the two cactus chloroplast genomes from NCBI RefSeq as guides, but so far I'm still waiting for the job to start.

### Annotate with PGA

I'm attempting this tool because GeSeq is hanging pretty badly...

```
git clone https://github.com/quxiaojian/PGA.git
PATH=$PATH:/data/biocore/programs/PGA
chmod -R u+rwx /data/biocore/programs/PGA
chmod -R g+rwx /data/biocore/programs/PGA

module load blast-plus-2.12.0-6z
perl PGA.pl -r test/angiosperms/reference -t test/angiosperms/target # successful test!
```

I loaded GenBank files for all Cactaceae chloroplast complete genomes available on NCBI, a total of 106 sequences. I'll use these as the reference files for PGA.

```
perl /data/biocore/programs/PGA/PGA.pl \
	-r /data/gencore/analysis_projects/8710018_Sunidhi/cactaceae-chloroplasts \
	-t /data/gencore/analysis_projects/8710018_Sunidhi/getOrganelle_output/pga
```

This output looks great - well-annotated GenBank formatted file. I'm going to create a bed file specifically for the Photosystem I genes (psaA, psaB, psaC, psaI, and psaJ - the others appear to be encoded in the nuclear genome, as do the Lhc subunits of the PSI-LHCI complex) so that I can extract their sequence from the fasta file.

BED files are considered to be 0-based (the first nucleotide in the chromosome is 0) while GenBank files are 1-based (the first nucleotide in the chromosome is 1). So, I'm using the guidance on this page: https://www.biostars.org/p/84686/ to convert accurately.

```
grep -B 1 'psa' embplant_pt.K105.complete.graph1.1.path_sequence_renamed.gb
     gene            complement(45762..48014)
                     /gene="psaA"
     CDS             complement(45762..48014)
                     /gene="psaA"
--
     gene            complement(43532..45736)
                     /gene="psaB"
     CDS             complement(43532..45736)
                     /gene="psaB"
--
     gene            118576..118821
                     /gene="psaC"
     CDS             118576..118821
                     /gene="psaC"
--
     gene            68883..68993
                     /gene="psaI"
     CDS             68883..68993
                     /gene="psaI"
--
     gene            complement(78780..78908)
                     /gene="psaJ"
     CDS             complement(78780..78908)
                     /gene="psaJ"
```

```
touch psa.bed
echo -e "chr1\t45761\t48014\tpsaA\t1\t-" >> psa.bed
echo -e "chr1\t43531\t45736\tpsaB\t1\t-" >> psa.bed
echo -e "chr1\t118575\t118821\tpsaC\t1\t+" >> psa.bed
echo -e "chr1\t68882\t68993\tpsaI\t1\t+" >> psa.bed
echo -e "chr1\t78779\t78908\tpsaJ\t1\t-" >> psa.bed

module load bedtools2-2.31.0-gw

bedtools getfasta -fi ../embplant_pt.K105.complete.graph1.1.path_sequence_renamed.fasta \
									-bed psa.bed -fo psa.fa

gcsplit psa.fa "/>/" "{4}" # then rename the individual files by their gene name
```

### Visualize Annotations with Chloroplot

The code for this is in the file chloroplot.R.

## Files to Return

### Input Data

#### Raw Fastq/QC

#### Cut/Trimmed Fastq/QC

#### getOrganelle Output (Params 1)

    * Optimized assembly files:
		    * graph1.1.path renamed fasta/fai
		    * graph1.1.path renamed gff3
				* graph1.selected gfa assembly graph
				* screen shots of graph and stats for the selected assembly graph

#### PGA Output (Annotation GenBank file)

#### Scripts

    * fastqc generation
		* cut-trim-filter
		* chloroplot.R
		* pga.sh
		* getOrganelle.sh
		* workflow overview

#### Write up methods/results overview with citations
