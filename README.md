# WGS_Assembly
Reference based and de novo assembly of whole genomes


# EANBiT_WGS_Assembly
Reference based and de novo assembly of whole genomes

## **Introduction**
This is a follow up on the series of bioinformatics trainings that happen at Uganda Virus Research Institute facilitated by H3ABioNet. This particular one focuses on the main procedure of metagenomics NGS data  analysis. Custom scripts are used but the users should take note that, the required software is available and accessible on respective machines. It may also be necessary to change paths to main directories to fit ones file system. This can always be changed at the top of each script.

## **Getting started**
Create a directory named "WGS" in your home directory. This is where the analysis will be carried out. As we go on, we shall be creating sub-directories accordingly to ensure that intermediate outputs can be refered to easily whenever needed for downstream analyses. 
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir WGS
```
##  **Reference based assembly**{#ref} 
Mapping to reference sequences is done using `bowtie2` [@bowtie2]. For this training, the references were indexed apriori and are available in the directory (`"$HOME/WGS/refs"`). First, make a directory `alignment` which will contain the alignment results. The script `align.sh` aligns reads to each of the reference sequences, converts .sam to .bam files and sorts them, index the sorted files, generates summary of mapped and unmapped reads and ultimately extracts mapped reads.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/alignment"
bowtie2-build -f $refname/$ref $refname/$refname
```

Map the reads onto the reference sequences
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
	bowtie2 -x $refdir/${refname}/${refname} -1 $fastqdir/${p}_R1*.fastq -2 $fastqdir/${p}_R2*.fastq -S $alignmnet/${p}_align/${p}_${refname}.sam >> $alignment/log.txt
```

Convert sam to bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
samtools view -bS $alignmnet/${p}_align/${p}_${refname}.sam > $alignmnet/${p}_align/${p}_${refname}.bam
```

Sort the bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
		samtools sort $alignmnet/${p}_align/${p}_${refname}.bam -o $alignmnet/${p}_align/${p}_${refname}.sort.bam
```
Only remain with the sorted bam file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
		mv $alignmnet/${p}_align/${p}_${refname}.sort.bam $alignmnet/${p}_align/${p}_${refname}.bam
```

Create an index for the alignment file
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
		samtools index $alignmnet/${p}_align/${p}_${refname}.bam
```

Get alignment stats of the alignment
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
        samtools idxstats $alignmnet/${p}_align/${p}_${refname}.bam
```

Obtain only mapped reads from the alignment file.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
        samtools view -bh -F4 $alignmnet/${p}_align/${p}_${refname}.bam > $alignmnet/${p}_align/${p}_${refname}.mapped.bam
```
##  **De novo assembly**{#denovo}
Here we use `spades` [@spades] to assemble host free reads to obtain consesus contigs and haplotype assembly.This script used paired end reads stored as fastq files. If the reads are not of this format, one may consider changing acordingly (by editing spades_assemble.sh). Among other outputs, the script produces consensus contigs and conservative regions of diploid genome.

```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/assembly"
fq2fa --merge  ${p}_L001_R1_001.fastq ${p}_L001_R2_001.fastq  ${p}_R12.fa
idba_ud -r ${p}_R12.fa -o ${p}
```
<!--
#### **Remove host sequences**
There is need to get rid of host DNA sequences that could have contaminated the data earlier in the sampling and extraction stage. This is done by mapping the reads to host reference genome and picking the unmapped sequences. At this point, we create a directory to store host-free sequence data. `Samtools` [@samtools] and `bedtools` [@bedtools] are the key tools used for this purpose.
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}
mkdir "$HOME/WGS/hostfree"

bowtie2 -x $hostgenome/host_DB -1 $trimmed/${p}_qtrim.1.fastq -2 $trimmed/${p}_qtrim.2.fastq --threads 20 -S $hostfree/${p}_mapped_and_unmapped.sam
        samtools view -bS $hostfree/${p}_mapped_and_unmapped.sam | samtools view -b -f 12 -F 256 | samtools sort -n > $hostfree/${p}_mapped_unmapped_sorted.bam
        rm -f $hostfree/${p}_mapped_and_unmapped.sam #remove sam files to save space
        bedtools bamtofastq -i $hostfree/${p}_mapped_unmapped_sorted.bam -fq ${p}_hostfree_R1.fastq -fq2 ${p}_hostfree_R2.fastq
        rm -f $hostfree/${p}_mapped_unmapped_sorted.bam 
```
-->

To evaluate and assess the assembly, we use `quast`. This will provide a summary of the metagenome assembly, including but not limited to N50, N75, L50, L75, GC percentage, number of contigs with size greater than 500bp (Only assesses the consensus, similar procedure can be used to assess other outputs).
```{r,eval=FALSE,error=FALSE,warning=FALSE,message=FALSE,echo=TRUE}

```
## **References**

