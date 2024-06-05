# Fish541-rnaseq
Project background: This project was part of the FISH 541 lab. The purpose of this project was to go through the process of RNA sequencing data analysis, starting with quality control, alignment (or psuedoalignment in this case), finding DEGs and finally functional analysis/gene ontology. 
The data used for this project was obtained from the Roberts lab nightingale directory. For the first half (until the DEG list step), atlantic salmon data were used. We then switched to seastar wasting disease data as it contained more samples (6 vs 2 samples) and allowed for better visualization for the purposes of this project. In my data and output directories, any YOJI file is related to atlantic salmon and the seastar samples are labled as "seastarwastingdisease". 
The Atlantic salmon samples were downloaded from https://owl.fish.washington.edu/nightingales/S_salar/. All the 8C_26psu samples were used for the first half of this project. 
For the seastar wasting disease samples, @Steven Roberts directly uploaded the abundance.tsv files of each sample into my Raven directory. These samples can be found in my data folder in this repository under "seastar abundance.tsv files"



The process was as follows: 
Step 1: Quality control
FastQC was run to assess the quality of the raw fastq file. A multiQC analysis was then be run on this output file to obtain a single consolidated report of all the fastQC files as fastQC performs the analysis on each individual sample. The output of multiQC is an html report file (multiqc_report.html).
```{bash}
/home/shared/FastQC-0.12.1/fastqc \
--threads 40 \
--outdir ../output/ \
/home/shared/8TB_HDD_02/davnoor/data/*.gz
```



Step 2: Alignment/pseudoalignment
Kallisto was used to perform pseudoaligmnent on salmo salar samples. A reference transcriptome was obtained from NCBI (rna.fna file that is uploaded onto my data folder) 
 This file should then be indexed to facilitate the k-mer generation that is required for Kallisto as follows:
```{r}
system(&quot;/home/shared/kallisto/kallisto index --index=&#39;kallisto_index_Salmosalar.idx&#39; ../data/rna.fna&quot;)
```
After indexing, align the sequence using the following code:
```{bash}
/home/shared/kallisto/kallisto quant \
-i ../code/kallisto_index_Salmosalar.idx \
-t 20 \
-o ../output/Kallistoquant_s_4_sequenceYOJI1 \
../data/s_4_sequenceYOJI1.txt.gz \
--single \
-l 50 \
-s 0.2
```
This step (alignment) was repeated for each of the (two) samples. The output was an
abundance.tsv file which should be saved in the output folder as instructed in the code. This file
contains abundance estimates (how much of a particular transcript is found in this sample, aka
expression level). Next, a trinity function was run to convert these abundance estimates into a
matrix for downstream processing as follows:
```{bash}
/home/shared/trinityrnaseq-v2.12.0/util/abundance_estimates_to_matrix.pl \
--est_method kallisto \
--gene_trans_map none \
--out_prefix ../output/Kallistoquant_s_4_sequenceYOJI1 \
--name_sample_by_basedir \
../output/Kallistoquant_s_4_sequenceYOJI1/abundance.tsv \
../output/Kallistoquant_s_4_sequenceYOJI2/abundance.tsv \
```
Step 3: Find DEGs (differentially expressed genes)
If the expression levels of transcripts from the abundance estimates found during
pseudoalignment is statistically different between the samples, the associated genes are
declared differentially expressed (4). To get a list of all these differentially expressed genes
(DEGs), DESeq2 was used in Rstudio. First, the DESeq2 package was installed using BiocManager as
follows:
```{r}
if (!require(&quot;BiocManager&quot;, quietly = TRUE))
install.packages(&quot;BiocManager&quot;)
BiocManager::install(&quot;DESeq2&quot;)
```
Load all the required functions from specified libraries as follows:
```{r}
library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(data.table)
```

Since a matrix of read counts was already prepared using trinity during the alignment step, the
count matrix was used as input for DESeq2. This code chunk processes the count matrix for
further analysis with DESeq2 by assigning row names, and removing the first column (which
was previously set to the row names:
```{r}
countmatrix &lt;- read.delim(&quot;../output/Kallistoquant_s_4_sequenceYOJI1.isoform.counts.matrix&quot;, header = TRUE,
sep = &#39;\t&#39;)
rownames(countmatrix) &lt;- countmatrix$X
countmatrix &lt;- countmatrix[,-1]
head(countmatrix)
```
The counts were rounded up to whole numbers:
```{r}
countmatrix &lt;- round(countmatrix, 0)
str(countmatrix)
```
A dataset was created for DESeq2 (note: I switched to seastar wasting disease analysis in this step):
```{r}
deseq2.colData &lt;- data.frame(condition=factor(c(rep(&quot;heatkilled&quot;, 1), rep(&quot;wastingdisease&quot;, 1))),
type=factor(rep(&quot;single-read&quot;, 2)))
rownames(deseq2.colData) &lt;- colnames(data)
deseq2.dds &lt;- DESeqDataSetFromMatrix(countData = countmatrix,
colData = deseq2.colData,
design = ~ condition)
```
Finally, the DESeq2 analysis was run:
```{r}
deseq2.dds &lt;- DESeq(deseq2.dds)
deseq2.res &lt;- results(deseq2.dds)
deseq2.res &lt;- deseq2.res[order(rownames(deseq2.res)), ]
```
And visualized using:
```{r}
head(deseq2.res)
```
The output of this analysis was a DEGlist.tab file in the output folder that contains a table
with a list of differentially expressed genes. This data was then used to plot various graphs such
as a PCA plot or MA plot. The code for these can be found in my code folder. 


Step 4: Functional analysis
Once a list of DEGs is obtained, it was used to obtain insights on what specific genes are
enriched and how these genes translate into specific processes and biological functions within
the organism. This was done with the help of Gene Ontology (GO), which is a database that
matches genes to their roles in three different categories: biological processes, molecular
function and cellular components. After generating GO terms, the DEGlist was matched to the GO
terms to find overrepresented categories to find which pathways are functionally enriched. To
do this, a UniProt Swiss-Prot (SP) database file was obtained from uniport.org that contains an annotated
protein sequence database using:
```{bash}
cd ../data
curl -O
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta.gz uniprot_sprot_r2024_02.fasta.gz
gunzip -k uniprot_sprot_r2024_02.fasta.gz
```
A directory was created using this database so that it is compatible with BLAST:
```{bash}
mkdir ../blastdb
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/uniprot_sprot_r2024_02.fasta \
-dbtype prot \
-out ../blastdb/uniprot_sprot_r2024_02
```
A BLASTX search was then run to identify protein sequences that corresponded to my
transcriptome. For the purpose of this lab, a transcriptome (Phel_transcriptome) was obtained from the
Roberts lab for seastar wasting disease and was output in a BLAST table format:
```{bash}
/home/shared/ncbi-blast-2.15.0+/bin/blastx \
-query ../data/Phel_transcriptome.fa \
-db ../blastdb/uniprot_sprot_r2024_02 \
-out ../output/Phel_uniprot_blastx.tab \
-evalue 1E-20 \
-num_threads 20 \
-max_target_seqs 1 \
-outfmt 6
```
Some reformatting of the table was required to match it with the DEGlist and translate it
into an R data frame (spgo), which was performed as follows:
```{bash}
tr &#39;|&#39; &#39;\t&#39; &lt; ../output/Phel_uniprot_blastx.tab \
&gt; ../output/Phel_uniprot_blastx_sep.tab
head -1 ../output/Phel_uniprot_blastx_sep.tab
```
```{r}
spgo &lt;- read.csv(&quot;https://gannet.fish.washington.edu/seashell/snaps/uniprot_table_r2023_01.tab&quot;, sep = &#39;\t&#39;,
header = TRUE)

install.packages(&#39;DT&#39;)
```
```{r}
blast &lt;- read.csv(&quot;../output/Phel_uniprot_blastx_sep.tab&quot;, sep = &#39;\t&#39;, header = FALSE)
```
Once the table was compatible with the R data frame, the columns were renamed since the
original table does not have column names. The first column name is formatted to be the same
as the DEG list first column name so that they can be merged later on:
```{r}
#Rename columns
cols &lt;- c(&quot;GeneID&quot;,&quot;SeqID&quot;,&quot;pident&quot;,&quot;length&quot;,&quot;mismatch&quot;,&quot;gapopen&quot;,&quot;qstart&quot;, &quot;qend&quot;, &quot;sstart&quot;, &quot;send&quot;, &quot;evalue&quot;,
&quot;bitscore&quot;)
colnames(blast) &lt;-cols
head(blast)
```
The necessary functions were loaded to merge the DEG list with the UniProt SP list. This is important as
the DEGs which only show the gene names will now be linked to GO terms which can give us
their functions. This will allowed me to obtain data about what specific processes/genes were
enriched and visualize this data. This was done as follows:
```{r}
library(tibble)
library(tidyverse)
```
```{r}
Joinedtable &lt;- merge(DEGlist, blast, by = &quot;GeneID&quot;)
head(Joinedtable)
```
The output of this was a joined table where the genes in the DEG list all are now associated with
their specific GO IDs. This list of GO IDs were then input into a functional annotation tool (DAVID bioinformatics) to obtain gene ontology information. In DAVID, the DEG list must be input as the gene list and the list of genes from the blast table must be used as the background. For this analysis UNIPROT_ACCESSION was used as the identifier. The output is a list of all the enriched genes sorted into the three GO categories: biological processes, cellular components and molecular functions. Each category has an output list of all the genes that were enriched for that particular category and their associated functions. Clicking on the output term leads to a QuickGo page that has the GO ID and a detailed flowchart of the processes that occur after the selected process (i.e., the downstream impacts of this particular gene being enriched) and their GO IDs. In this project, to find what gene the GO ID corresponds to, a blastquery GOslim table made by Grace Crandall was used. This table has a list of all the relevant genes matched with their corresponding GO IDs. The blastquery table was merged it with my DEG list to obtain all the GO IDs that match my gene list as follows:
```{r}
Blastquery &lt;- read.csv(&quot;../output/blastquery_GOslim/Blastquery-GOslim.txt&quot;, sep = &#39;\t&#39;, header = FALSE)
```
```{r}
#Rename columns
cols &lt;- c(&quot;GeneID&quot;,&quot;GO_ID&quot;,&quot;biological_process&quot;,&quot;category&quot;)
colnames(Blastquery) &lt;-cols
head(Blastquery)
```

```{r}
Joinedtable2 &lt;- merge(DEGlist, Blastquery, by = &quot;GeneID&quot;)
head(Joinedtable2)
```
```{r}
write.table(Joinedtable2, &quot;../output/Joinedtable2.tab&quot;, sep = &#39;\t&#39;, row.names = T)
```
This last step made it possible to obtain the names of the genes associated with the GO IDs.
The genes can now be looked up in Uniprot or a similar software to get more information about
the pathways that were enriched in the experimental sample as compared to the control.
