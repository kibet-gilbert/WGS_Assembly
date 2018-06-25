## **Introduction**
This is a follow up on the series of bioinformatics training sessions for the EAST African Network for Bioinformatics Trainin (EANBiT) held at KEMRI wellcome Trust, Kilifi Campus. This session focuses on the main procedure of  WGS reference based and de novo assembly.

NOTE: It may  be necessary to make slight changes along the way.

## **Getting started**
Create a directory named "WGS" in your home directory and navigate to it. This is where the analysis will be carried out. As we go on, we shall be creating sub-directories accordingly to ensure that intermediate outputs can be refered to easily whenever needed for downstream analyses. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir WGS
```

Inside WGS, make a directory `data` to store data. Download the data from a repository using `wget` program. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/data"
cd data
wget https://github.com/AlfredUg/WGS_Assembly/raw/master/reads_1.fq
wget https://github.com/AlfredUg/WGS_Assembly/raw/master/reads_2.fq
```

To confirm that the links have been successfully created, list the contents of the folder using the `ls` command. If the data is compressed (i.e .gz files), uncompress them using `gunzip *`

## **Quality control check**
Here we use `fastqc`, it is necessary that we store quality control files for easy reference. In WGS, create a sub-directory qcresults, this is where the fastqc results will be stored.
After that, do the quality checks;
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/qcresults"
fastqc "$HOME/WGS/data/*" -o "$HOME/WGS/qcresults"
```
Once this is done, navigate to `qcresults` and download the ".html" files to local machine and open them in any browser. This report can be used to assess quality distribution, length distribution, GC-content, nucleotide distribution. This informs downstream analysis.

## **Data quality trimming**
We use `trim_galore` for quality and adapter trimming. Depending on the qc results, it would be necessary to change some parameters used here accordingly. This script clips off the first 16 bases of the reads from the 5' end. In addition, it removes bases with phred quality less than 25 on the 3' end of the reads. We need to store quality trimmed reads, as such, in WGS directory, make a sub directory `trimmed`.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/trimmed"
trim_galore -q 25 -l 75 --dont_gzip --clip_R1 16 --clip_R2 16 --paired read_R1.fq read_R2.fq -o trimmed
```
If all goes well, trimmed reads will be available in trimmed.
You may consider looking at the trimmed reads using `fastqc` to check the improvement made by `trim_galore`.

##  **Reference based assembly** 

Now that we have quality reads, we can proceed to map the reads onto the reference genome. Here, we use reference `bowtie2`. Make a directory `genome` to store the reference genome.  Create a soft-link of the genome to this directory.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/reference"
cd reference
wget https://github.com/AlfredUg/WGS_Assembly/raw/master/lambda_virus.fa
```

Create a directory `alignment` that will contain assembly results. Here we use the mapping tool known as bwa. Navigate to alignement and perform the assembly. This proceeds in two steps, first is to create align  the reads on the reference and then summarising the alignment in a SAM file.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir alignment
cd alignment
bwa aln ../reference/lambda_virus.fa ../data/reads_1.fq > reads_1.sai
bwa aln ../reference/lambda_virus.fa ../data/reads_2.fq > reads_2.sai
bwa sampe ../reference/lambda_virus.fa reads_1.sai reads_2.sai ../data/reads_1.fq ../data/reads_2.fq > reads12_alignment.sam
```

Use command `ls` to view the contents of reference directory.

Convert sam to bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bS reads12_alignment.sam > reads12_alignment.bam
```

Sort the bam file such that it can easily be indexed.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools sort reads12_alignment.bam -o reads12_alignment_sorted.bam
```
Note: At this point, one can discard the SAM file, since we already have a lighter version of it. If need be, we can always convert back to a SAM file.

Create an index and compute stats for the alignment
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools index reads12_alignment_sorted.bam
samtools idxstats reads12_alignment_sorted.bam
```

One can choose to obtain only the mapped or unmapped reads from the entire alignment file.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bh -F4 reads12_alignment.bam > reads12_alignment.bam
```

Download the alignment file and open it in your favorite editor. See what each of the columns contain and mean.

Visualise the alignment using IGV.

##  **De novo assembly**
Here we use `velvet` to assemble reads to obtain consesus contigs. This proceedes in two stages, first is to create the hash index (using velveth) and then the assembly (velvetg).
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/denovo"
velveth denovo 29 read1.fa read2.fa -short
velvetg denovo
```

To evaluate and assess the assembly, we use `quast`. This will provide a summary of the metagenome assembly, including but not limited to N50, N75, L50, L75, GC percentage, number of contigs with size greater than 500bp (Only assesses the consensus).
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
./quast.py contigs_1.fasta  -R test_data/reference.fasta.gz \
        -G test_data/genes.txt \
        -O test_data/operons.txt \
        -1 test_data/reads1.fastq.gz -2 test_data/reads2.fastq.gz \
        -o quast_test_output
```
## **References**
