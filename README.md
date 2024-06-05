# Fish541-rnaseq
#Note: All code for this experiment can be found under the code folder 
#Step 1: Quality control
#FastQC was run to assess the quality of the raw fastq file. A multiQC analysis was then be run on this output file to obtain a single consolidated report of all the fastQC #files as fastQC performs the analysis on each individual sample. The output of multiQC is an html report file (multiqc_report.html).
#Step 2: Alignment/pseudoalignment
#Kallisto was used to perform pseudoaligmnent on salmo salar samples. A reference transcriptome was obtained from NCBI (rna.fna file that is uploaded onto my data folder) 
 This file should then be indexed to facilitate the k-mer generation that is
required for Kallisto as follows:
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
This step (alignment) should be repeated for each of the samples. The output will be an
abundance.tsv file which should be saved in the output folder as instructed in the code. This file
contains abundance estimates (how much of a particular transcript is found in this sample, aka
expression level). Next, a trinity function is run to convert these abundance estimates into a
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
(DEGs), DESeq2 will be used in Rstudio. First, install the DESeq2 package using BiocManager as
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
count matrix can be used as input for DESeq2. This code chunk processes the count matrix for
further analysis with DESeq2 by assigning row names, and removing the first column (which
was previously set to the row names:
```{r}
countmatrix &lt;- read.delim(&quot;../output/Kallistoquant_s_4_sequenceYOJI1.isoform.counts.matrix&quot;, header = TRUE,
sep = &#39;\t&#39;)
rownames(countmatrix) &lt;- countmatrix$X
countmatrix &lt;- countmatrix[,-1]
head(countmatrix)
```
Round the counts up to whole numbers:
```{r}
countmatrix &lt;- round(countmatrix, 0)
str(countmatrix)
```
Create a dataset for DESeq2 (note: I switched to seastar wasting disease analysis in this step):
```{r}
deseq2.colData &lt;- data.frame(condition=factor(c(rep(&quot;heatkilled&quot;, 1), rep(&quot;wastingdisease&quot;, 1))),
type=factor(rep(&quot;single-read&quot;, 2)))
rownames(deseq2.colData) &lt;- colnames(data)
deseq2.dds &lt;- DESeqDataSetFromMatrix(countData = countmatrix,
colData = deseq2.colData,
design = ~ condition)
```
Finally, run the DESeq2 analysis:
```{r}
deseq2.dds &lt;- DESeq(deseq2.dds)
deseq2.res &lt;- results(deseq2.dds)
deseq2.res &lt;- deseq2.res[order(rownames(deseq2.res)), ]
```
And visualize using:
```{r}
head(deseq2.res)
```
The output of this analysis will be a DEGlist.tab file in the output folder that contains a table
with a list of differentially expressed genes. This data can then used to plot various graphs such
as a PCA plot or MA plot as shown below. The code for these can be found on my Github
repository (davnoor/Fish541-rnaseq).

Figure 1: A: MA plot of DEGs indicating genes upregulated vs downregulated in seastars
affected by desiccation. B: PCA plot of genes expressed by control (heatkilled) group of seastars
vs seastars affected by wasting disease.
Step 4: Functional analysis
Once a list of DEGs is obtained, it can be used to obtain insights on what specific genes are
enriched and how these genes translate into specific processes and biological functions within
the organism. This can be done with the help of Gene Ontology (GO), which is a database that
matches genes to their roles in three different categories: biological processes, molecular
function and cellular components. After generating GO terms, match the DEGlist to the GO
terms to find overrepresented categories to find which pathways are functionally enriched. To
do this, use a UniProt Swiss-Prot (SP) database file from uniport.org that contains an annotated
protein sequence database using:
```{bash}
cd ../data
curl -O
https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
mv uniprot_sprot.fasta.gz uniprot_sprot_r2024_02.fasta.gz
gunzip -k uniprot_sprot_r2024_02.fasta.gz
```
Create a directory using this database so that it is compatible with BLAST:
```{bash}
mkdir ../blastdb
/home/shared/ncbi-blast-2.11.0+/bin/makeblastdb \
-in ../data/uniprot_sprot_r2024_02.fasta \
-dbtype prot \
-out ../blastdb/uniprot_sprot_r2024_02
```
A BLASTX search should then be run to identify protein sequences that correspond to your
transcriptome. In this example, a transcriptome (Phel_transcriptome) was obtained from the
Roberts lab for seastar wasting disease and output it in a BLAST table format:
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
Some reformatting of the table might be required to match it with the DEGlist and translate it
into an R data frame (spgo), which can be performed as follows:
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
Once the table is compatible with the R data frame, the columns must be renamed since the
original table does not have column names. The first column name is formatted to be the same
as the DEG list first column name so that they can be merged later on:
```{r}
#Rename columns
cols &lt;- c(&quot;GeneID&quot;,&quot;SeqID&quot;,&quot;pident&quot;,&quot;length&quot;,&quot;mismatch&quot;,&quot;gapopen&quot;,&quot;qstart&quot;, &quot;qend&quot;, &quot;sstart&quot;, &quot;send&quot;, &quot;evalue&quot;,
&quot;bitscore&quot;)
colnames(blast) &lt;-cols
head(blast)
```
Load the necessary functions to merge the DEG list with the UniProt SP list. This is important as
the DEGs which only show the gene names will now be linked to GO terms which can give us
their functions. This will allow us to obtain data about what specific processes/genes were
enriched and visualize this data. This can be done as follows:
```{r}
library(tibble)
library(tidyverse)
```
```{r}
Joinedtable &lt;- merge(DEGlist, blast, by = &quot;GeneID&quot;)
head(Joinedtable)
```
The output of this is a joined table where the genes in the DEG list all are now associated with
their specific GO IDs. This list of GO IDs can now be input into a functional annotation tool like
DAVID bioinformatics to obtain gene ontology information. In DAVID, the DEG list must be input
as the gene list and the list of genes from the blast table must be used as the background. For
this analysis UNIPROT_ACCESSION was used as the identifier. The output is a list of all the
enriched genes sorted into the three GO categories: biological processes, cellular components
and molecular functions. Each category has an output list of all the genes that were enriched
for that particular category and their associated functions. Clicking on the output term leads to
a QuickGo page that has the GO ID and a detailed flowchart of the processes that occur after
the selected process (i.e., the downstream impacts of this particular gene being enriched) and
their GO IDs. In this example, to find what gene the GO ID corresponds to, a blastquery GOslim
table made by Grace Crandall was used. This table has a list of all the relevant genes matched
with their corresponding GO IDs. The blastquery table was merged it with my DEG list to obtain
all the GO IDs that match my gene list as follows:
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
This last step makes it possible to obtain the names of the genes associated with the GO IDs.
The genes can now be looked up in Uniprot or a similar software to get more information about
the pathways that were enriched in the experimental sample as compared to the control.
