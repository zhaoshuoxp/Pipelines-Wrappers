# ATACseq, ChIPseq, RNAseq pipeline wrappers and more

-----
This repository has the following combined shell/awk/python/R scripts which can be used for High-throughput sequecning data analysis.

 * [ATACseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#atacseqsh): bulk ATACseq pipeline, from fastq to open chromatin regions.
 * [ChIPseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#chipseqsh): ChIPseq pipeline, from fastq to peak calling step.
 * [RNAseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rnaseqsh): bulk RNAseq pipeline, from fastq to differential expression genes.
 * [adapt_trim.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#adapt_trimsh): adapter trimming function, seperated from the above pipelines.
 * [cisVar.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#cisvarsh): pipeline wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).
 * [GATK_HF.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#gatk_hfsh): variants calling by [GATK](https://software.broadinstitute.org/gatk/), from fastq to vcf.
 * [trans_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#trans_assemblesh): *de novo* transcript assembly, from fatsq to GTF.
 * [PLAR.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#plarsh): *de novo* [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) lncRNA discovery pipeline wrapper.
 * [rRNA_dep.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rrna_depsh): ribosomal RNA depletion from fastq files.

> Requirements:
> Python3, cutadapt, macs2(>=2.1.1), R, DESeq2, featureCounts, bowtie2, bwa,STAR, fastqc, samtools, bedtools, bedGraphToBigWig, bedItemOverlapCount



-----

## ATACseq.sh

This script QC fastq files and align reads to hg19/GRCh37 using Bowtie2, convert to filtered BAM/BED and bigwig format, then call peaks with MACS2 in BEDPE mode after Tn5 shifting. 

#### Input
Paired-end fastq files with **_R1/2** extension, ex:test_R1.fastq.gz, test_R2.fastq.gz 
> Single-end sequencing data is also supported with -s, although it is not recommended.

#### Options

help message can be shown by `ATACseq.sh -h`

```shell
Usage: ATAC.sh <options> <reads1>|..<reads2> 
    Options:
        -i [str] Bowtie2 index PATH
        -p [str] Prefix of output
        -t [int] Threads (1 default)
        -s Single-end mod (DO NOT recommend, Paired-end default)
        -h Print this help message
```

#### Example

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ATACseq.sh
chmod 755 ATACseq.sh
./ATACseq.sh -i /path/to/bwt2idx/ -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

####  Output
All results will be store in current (./) directory.

* test_trimmed_R1/2.fastq.gz: adapter trimmed fastq files.

* test_mkdup.bam: all alignments, with duplicates marked.

* test_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.

* test_se.bed: useful filtered alignments in BED format.

* test_pe.bed: useful filtered alignments in BEDPE format, the 2nd and 3rd columns indicate the fragment start and end coordinates on genome.

* test.bw: bigwig file converted from test_se.bed, can be upload to genome browser for visualization.

* test_shift.bed: Tn5 shifted BEDPE format, it will be used for macs2 peak calling.

* macs2: output of macs2, see [here](https://github.com/taoliu/MACS#output-files). Only broad peaks will bed called by defualt. In addtion, test_broad_filtered.bed is the peaks file with hg19 blacklist filtering.

* fastqc: the report(s) of fastqc

* logs: running logs

  

-----

## ChIPseq.sh

This script QC fastq files and align reads to hg19/GRCh37(depends on index and GTF provided) using BWA, convert to filtered BAM/BED and bigwig format but DOES NOT call peaks.

#### Input
Paired-end fastq files with **_R1/2** extension, ex:test_R1.fastq.gz, test_R2.fastq.gz 
Or single-end fastq file with -p.

#### Options
help message can be shown by `ChIPseq.sh -h`

```shell
Usage: ChIPseq.sh <options> <reads1>|..<reads2> 
    Options:
        -i [str] BWA index PATH
        -p [str] Prefix of output
        -t [int] Threads (1 default)
        -s Single-end mod (Paired-end default)
        -a Use BWA aln algorithm (BWA mem default)
        -h Print this help message
```

#### Example

```shell
    wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ChIPseq.sh
    chmod 755 ChIPseq.sh
    ./ChIPseq.sh -i /path/to/bwaidx/ -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```
####  Output
All results will be store in current (./) directory.

* test_trimmed[_R1/2].fastq.gz: adapter trimmed fastq files.
* test_mkdup.bam: all alignments, with duplicates marked.
* test_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.
* test_se.bed: useful filtered alignments in BED format.
* test_pe.bed: useful filtered alignments in BEDPE format, the 2nd and 3rd columns indicate the fragment start and end coordinates on genome.
* test.bw: bigwig file converted from test_se.bed, can be upload to genome browser for visualization.
* fastqc: the report(s) of fastqc
* logs: running logs


#### Peak calling
> NOTE:
> this pipeline does NOT call peaks, you might want to run it manually.
> input is highly recommended for peak calling, put input fastq files through this pipeline with same parameter(s).

test_pe.bed (and input_pe.bed) can be used for macs2 peak calling in BEDPE mode:

```shell
macs2 callpeaks -t test_pe.bed -c input_pe.bed -f BEDPE -g hs -n test -B --SPMR
```

> --broad is recommended for histone modifications when using macs2

test_se.bed and test_filtered.bam can also be used in BED or BAM mode of macs2.

See more about [MACS2](https://github.com/taoliu/MACS) (for TFs peak calling) and [SICER](https://home.gwu.edu/~wpeng/Software.htm) or [SICERpy](https://github.com/dariober/SICERpy) (for Histone Mods peak calling).



-----

## RNAseq.sh

This script QC fastq files and align reads to hg19/GRCh37(depends on index and GTF provided) using STAR, featureCounts and DESeq2 will be used for reads count and differntial expresssion genes discovery.

#### Input
* Paired-end fastq files with **_R1/2.fastq.gz** extension, put fastq files of each condition all together in a directoy, e.g.

> Single-end not supported

```shell
ls -1 ./
cond1_rep1_R1.fastq.gz 
cond1_rep1_R2.fastq.gz 
cond1_rep2_R1.fastq.gz 
cond1_rep2_R2.fastq.gz
cond2_rep1_R1.fastq.gz 
cond2_rep1_R2.fastq.gz 
cond2_rep2_R1.fastq.gz 
cond2_rep2_R2.fastq.gz
....
```

And a text file discribing samples per conditon e.g.

```shell
sample  condition
cond1_rep1  cond1
cond1_rep2  cond1
cond2_rep1  cond2
cond2_rep2  cond2
....
```

You can use the script to scan fastq files and generate the text file:

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/RNAseq.sh
chmod 755 RNAseq.sh 
./RNAseq.sh -p /path/to/directory/contains/fastq/
```
then the condition.txt will be created and open with VIM. sample column will have been filled, edit the text by adding the condition information on 2nd column.

> NOTE:
> !!! Provide **the PATH of the DIRECTORY** which contains fastq to the scripts, DO NOT give the path of fastq files directly !!!

#### Options
help message can be shown by `RNAseq.sh -h`

```shell
Usage: RNAseq.sh <options> -c conditions.txt /PATH/to/directoy/contains/fastq/ 
    Options:
        -c [str] /PATH/to/conditions.txt
        -i [str] STAR index PATH
        -g [str] Reference GTF transcripts PATH
        -t [int] Threads (1 default)
        -p prepare condition.txt
        -n Nextera adapters (Truseq default)
        -h Print this help message
```

#### Example

```shell
./RNAseq.sh -i /path/to/STARidx/ -g /path/to/ref/GTF -c conditions.txt -t 24 /path/to/directory/contains/fastq/
```

####  Output

All results will be store in current (./) directory.

* {prefix}_R1/2_trimmed.gz: adapter trimmed fastq files.
* {prefix}.bam: STAR output, accepted alignments.
* count.txt: featureCounts output, raw fragments count.
* deseq.r: R script for DESeq2.
* all_genes_exp.txt: size-factor normalized gene expression levels with *P* values.
* fastqc: the report(s) of fastqc
* logs: running logs

> NOTE:
 Sample names in conditions.txt has to match featureCount output, check your text or generate it by the script.
 This script cannot compare the DE genes condition by condition automatically if you have >2 conditions to compare. Either edit deseq.r or load count.txt to R. A online tool can be used [iDEP](http://bioinformatics.sdstate.edu/idep/).



----

## adapt_trim.sh

This script is seperated from ChIPseq.sh, it trims adapter sequences from fastq files with cutadapt@python3.

#### Input
* Paired-end fastq files or single-end with -s, e.g. test_R1.fastq.gz test_R2.fastq.gz

#### Options
help message can be shown by `adapt_trim.sh -h`

```shell
Usage: adapt_trim.sh <options> <reads1>|..<reads2
    Options:
    	-p Prefix of output
        -t Threads (1 default)
        -s Single-end mod (Paired-end default)
        -n Nextera adapters (Truseq default)
        -h Print this help message
```

#### Example

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/adapt_trim.sh
chmod 755 adapt_trim.sh 
./adapt_trim.sh -p test -t 24 test_R1.fastq.gz test_R2.fastq.gz
```

> NOTE:
multi-threads support only works with python3>=3.4, multiprocessing>=0.70, cutadapt>=1.15 and pigz

#### Output
All results will be store in current (./) directory.

* test_R1/2_trimmed.gz: adapter trimmed fastq files.

  

------
## trans_assemble.sh

This script QC fastq files and align reads to hg19/GRCh37(depends on index and GTF provided) using HISAT2. *De novo* transcripts assembly will be performed by stingtie.

#### Input
* Paired-end fastq files.

#### Usage

```shell
./trans_assemble.sh <reads1> <reads2> <prefix of output> <starnd: fr|rf|un>
```

> NOTE:
Edit the script and mofiy $threads, $index, $gtf.

#### Output

All results will be store in current (./) directory.

* {prefix}.bam: sorted accepted alignments.

* {prefix}.gtf: *de novo* transcripts assembled with reference GTF guiding.

* logs: running logs.

  

------
## cisVar.sh

This script is a wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).

#### Input
* Output of [hornet](https://github.com/TheFraserLab/Hornet), sorted and indexed bam file.

#### Usage

```shell
./cisVar.sh hornet.bam <read_depth> <indivdual files>
```

> NOTE:
Edit the script and mofiy $vcf PATH to your SNP vcf files.

See more about [cisVar](https://github.com/TheFraserLab/cisVar).

#### Output
All results will be store in current (./) directory.

* {prefix}.{read_depth}.final.txt: mian regression output.

* {prefix}.${DEP}.total.txt.prepost.png: desity plot of the regression output.

  

-----

## PLAR.sh

This script is a re-wrrited wrapper of [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR).  
> Requirements
stringtie, cuffdiff, plar, CPC2 and HMMER

#### Input
* BAM files of aligned RNAseq data, and coresponding GTFs.

#### Usage

```shell
./PLAR.sh hornet.bam <output dir> <prefix_sample1,prefix_sample2...> <strand:rf|fr|un> <sample1_rep1.gtf> <sample1_rep2.gtf> <sample2_rep1.gtf> ... <sample1_rep1.bam,sample1_rep2.bam,sample2_rep1..>
```

> NOTE:
Edit the script and mofiy $plar_path, $cpc2_path.
Addtional annotation files required in $plar_path, see [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) for nore details.

#### Output
All results will be store in current (./) directory.

* final_lncRNA.bed: filtered and sorted lncRNAs in BED12 format.

  

------
## rRNA_dep.sh

This script removes ribosomal RNA reads from fastq files by mapping them to rRNA genes and retrieving unmmaped reads.

#### Input
* Paired-end fastq files.

#### Usage

```shell
./rRNA_dep.sh <reads1> <reads2> <prefix of output>
```

> NOTE:
Edit the script and mofiy $threads, $genome, $gtf.

#### Output

All results will be store in current (./) directory.

* {prefix}_dep_R1/2_fq.gz: rRNA removed fastq files.

* {prefix}_rRNA.log: mapping log.

  

------

## GATK_HF.sh

This script is a wrapper of variants calling by [GATK](https://software.broadinstitute.org/gatk/).

#### Input
* Paired-end fastq files.

#### Usage

```shell
./GATK_HF.sh fastq1 fastq2
```

> NOTE:
Edit the script and mofiy $hg19, $picard, $gatk, $gatk_bundle_hg19, $gatk_ref_hg19 PATH.

See more about [GATK Best Practices](https://software.broadinstitute.org/gatk/best-practices/).

#### Output

All results will be store in current (./) directory.

* GVCF and VCF

  

------
Author [@zhaoshuoxp](https://github.com/zhaoshuoxp)  
Mar 27 2019  