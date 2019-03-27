# ATACseq, ChIPseq, RNAseq pipeline wrappers amd more
-----
This repository has the following combined shell/awk/python/R scripts which can be used for High-throughput sequecning data analysis.

 * ATACseq.sh: bulk ATACseq pipeline, from fastq to open chromatin regions.
 * ChIPseq.sh: ChIPseq pipeline, from fastq to peak calling step.
 * RNAseq.sh: bulk ATACseq pipeline, from fastq to differential expression genes.
 * adapt_trim.sh: adapter trimming function, seperated from the above pipelines.
 * cisVar.sh: pipeline wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).
 * GATK_HF.sh: variants calling by [GATK](https://software.broadinstitute.org/gatk/), from fastq to vcf.
 * trans_assemble.sh: *de novo* transcript assembly, from fatsq to GTF.
 * PLAR.sh: *de novo* [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) lncRNA discovery pipeline wrapper.
 * rRNA_dep.sh: ribosomal RNA depletion from fastq files.

> Requirements:
> Python3, cutadapt, macs2(>=2.1.1), R, DESeq2, featureCounts, bowtie2, bwa,STAR, fastqc, samtools, bedtools, bedGraphToBigWig, bedItemOverlapCount

-----

### ATACseq.sh

This script will QC fastq files and align to hg19/GRCh37 using Bowtie2, convert to filtered BAM/BED and bigwig format, then call peaks with MACS2 in BEDPE mode after Tn5 shifting. All results will be store in current (./) directory.

#### Input
Paired-end fastq files with _R1/2 extension, ex:test_R1.fastq.gz, test_R2.fastq.gz 
> Single-end sequencing data is also supported with -s, although it is not recommended.

#### Options

help message can be shown by `ATACseq.sh -h`

    Usage: ATAC.sh <options> <reads1>|..<reads2> 
        -i [str] Bowtie2 index PATH
        -p [str] Prefix of output
        -t [int] Threads (1 default)
        -s Single-end mod (DO NOT recommend, Paired-end default)
        -h Print this help message

#### Example run

    wget https://github.com/zhaoshuoxp/Pipelines-Wrappers/blob/master/ATACseq.sh
    chmod 755 ATACseq.sh
    ./ATACseq.sh -i /path/to/bwt2idx/ -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz

####  Output

* ${prefix}_trimmed_R1/2.fastq.gz: adapter trimmed fastq files.
* ${prefix}_mkdup.bam: all alignments, with duplicates marked.
* ${prefix}_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.
* 
