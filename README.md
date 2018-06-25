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
Here we use `fastqc` [@fastqc], it is necessary that we store quality control files for easy reference. In WGS, create a sub-directory qcresults, this is where the fastqc results will be stored.
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
trim_galore -q 25 -l 75 --dont_gzip --clip_R1 16 --clip_R2 16 --paired read_R1.fastq read_R2.fastq -o trimmed
```
If all goes well, trimmed reads will be available in trimmed.
Below is the details of trim.sh script, it iteratively trims all fastq files in the `data` directory. You may consider looking at the trimmed reads using `fastqc` to check the improvement made by `trim_galore`.

##  **Reference based assembly** 

Now that we have quality reads, we can proceed to map the reads onto the reference genome. Here, we use reference `bowtie2`. Make a directory `genome` to store the reference genome.  Create a soft-link of the genome to this directory.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/reference"
ln -s "path/to/reference.fa" "/data"
```

Create a directory that will contain assembly results.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "HOME/WGS/alignment" 
bowtie2-build -f "path/to/reference" "path/to/reference/basename"
```

Map the reads onto the reference sequences. Look at the documentation of bowtie2 to understand and see if you make any changes.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
bowtie2 -x path/to/reference -1 path/to/read1.fq -2 path/to/read2.fq -S path/to/alignmnet/alignment.sam >> path/to/alignment/log.txt
```

Convert sam to bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bS alignment.sam > alignmnet.bam
```

Sort the bam file such that it can easily be indexed.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools sort alignment.bam -o alignmnet_sort.bam
```
Note: At this point, one can discard the sam file, since we already have a lighter version of it. If need be we can always convert back to a sam file.

Create an index and compute stats for the alignment
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools index alignmnet.bam
samtools idxstats alignmnet.bam
```

One can choose to obtain only the mapped or unmapped reads from the entire alignment file.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bh -F4 alignmnet.bam > $alignmnent_mapped.bam
```

Download thealignment file and open it in your favorite editor. See what each of the columns contain and mean.

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

```
## **References**
