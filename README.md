# Single cell and bulk sequencing pipeline wrappers

-----
This repository has the following combined shell/awk/python/R scripts which can be used for High-throughput sequecning data analysis.

* [cellranger.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#cellrangersh): Run CellRanger in batches to process 10x Genomics single-cell fastq files.
 * [ATACseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#atacseqsh): bulk ATACseq and CUT&TAG pipeline, from fastq to open chromatin regions/peaks.
 * [ChIPseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#chipseqsh): ChIPseq pipeline, from fastq to peak calling step.
 * [RNAseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rnaseqsh): bulk RNAseq pipeline, from fastq to differentially expressed genes.
 * [adapt_trim.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#adapt_trimsh): adapter trimming function, seperated from the above pipelines.
 * [cisVar.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#cisvarsh): pipeline wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).
 * [trans_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#trans_assemblesh): *de novo* transcript assembly, from fastq to GTF.
 * [PLAR.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#plarsh): *de novo* [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) lncRNA discovery pipeline wrapper.
 * [rRNA_dep.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rrna_depsh): ribosomal RNA depletion from fastq files.
 * [CRISPRlib.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#crisprlibsh): mapping CRISPR sgRNA library, from fastq to tables.

> Requirements:
> Python3, cutadapt, macs2(>=2.1.1), R, DESeq2, featureCounts, bowtie1/2, bwa,STAR, fastqc, samtools, bedtools, deeptools, cellranger(7 and 8), cellranger-atac and cellranger-arc

[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu) [![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

-----

## cellranger.sh

This script runs cellranger automatically for all samples in the specified fastq directory. Sample names are inferred from the fastq filenames, i.e`<prefix>_S*_L*_R?_00*.f*q.gz`. Outputs are saved in the current directory with filenames based on the inferred prefix. The script supports 10x genomic RNA, ATAC, and multiome based on the specified data type (-m). The script can also perform aggregations for ATAC and multiome data using `-a` , and optionally `-c` to designate an aggregation csv file. You can find more details on setting up the aggregation csv file for [ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/aggr#csv_setup) and [multiome](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/aggregating-multiple-gem-wells-aggr#setting-up-an-aggregation-csv).

> Data type (-m) is also required for aggregation (-a). If no aggregation csv file designated, a csv file will be automatically generated and ALL folders/samples in the provided path will be aggregated.

#### Input

Paired-end FASTQ files following the cellranger demultiplexed naming conventions (cellranger mkfastq), all samples together in a directory.

#### Options

help message can be shown by `cellranger.sh -h`

```shell
 Usage: cellranger.sh [options] <fastq_directory | output_directory>

Options:
  -g [str]    Genome build <hg38|mm10> (required unless -x is set)
  -m [str]    Data type <rna|atac|multiome> (required)
  -x [str]    Custom reference path (overrides -g)
  -t [int]    Threads (default: 20)
  -r [int]    Memory in GB (default: 200)
  -u [int]    Use CellRanger version <8|7> (for -m rna only)
  -a          Run aggregation (aggr) mode
  -c [str]    Custom CSV file for aggr (optional)
  -n          Normalize in aggr mode (default: none)
  -s          Enable secondary analysis (default: off)
  --gex_path  RNA fastq path (for multiome mode)
  --atac_path ATAC fastq path (for multiome mode)
  -h          Show help
```

#### Example

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/cellranger.sh
chmod 755 cellranger.sh
# run cellranger count individually
./cellranger.sh -g mm10 -m atac /path/to/all_samples/
# run cellranger aggr
./cellranger.sh -g mm10 -m atac -a /path/to/cellranger/count/output/with/all/samples/

# run cellranger multiome individually
./cellranger.sh -g mm10 -m multiome --gex_path /path/to/all_samples/GEX/ --atac_path /path/to/all_samples/ATAC/
```

> Sequencing depth normalization is not recommended in most situations. For RNA data, aggregating the data matrix in Seurat/Scanpy is a more efficient way.

####  Output

All results will be store in current (./) directory.

* For indivdual cellranger count run, each samples will be stored in a filder named by their predicted prefixes.
* For cellranger aggr run, aggregated data will be stored in "aggr" folder and a aggr.csv will be created if no csv file designated.

> For more details of the cellranger output, see [RNA](https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/outputs/cr-outputs-overview), [ATAC](https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/using/count), and [Multiome](https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/analysis/outputs/understanding-output)



-----

## ATACseq.sh

This script checks quality controls of fastq files, then aligns reads to the specified reference genome using Bowtie2 or chromap, depending on the selected species passed by -g or the provided index and other necessary files specified by -i, -b, and -c. It converts the alignments to filtered BAM/BED and bigwig formats, and subsequently identifies peaks using MACS2 in BEDPE mode following Tn5 shifting.

> This script works for both ATACseq and CUT&TAG.

#### Input
Paired-end fastq files with **_R1/2** suffix, i.e. test_R1.fastq.gz, test_R2.fastq.gz 
> Single-end sequencing data is also supported with -s, but it is not recommended.

#### Options

help message can be shown by `ATACseq.sh -h`

```shell
Usage: ATAC.sh <options> -g <hg38|hg19|mm10> <reads1>|..<reads2> 
  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -x [str] Custom bowtie2 index PATH
    -b [str] Custom blacklist PATH
    -m [str] Genome size abbr supported by MACS2
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -s Single-end mod (DO NOT recommend, Paired-end default)
    -c Using chromap to process FASTQ instead of canonical bowtie2
    -i [str] Custom chromap genome index (only valid with -c option)
    -r [str] Custom chromap genome reference (only valid with -c option)
    -z [str] Custom chromosome size table
    -h Print this help message
```

#### Example

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ATACseq.sh
chmod 755 ATACseq.sh
./ATACseq.sh -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

Alternatively, you may use [chromap](https://github.com/haowenz/chromap) aligner to speed up the processing:

```shell
./ATACseq.sh -c -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

> chromap is recomended for most cases, as it is ultra fast and outputs a little bit less (<5%) aligned reads.

####  Output

All results will be store in current (./) directory.

* {prefix}_trimmed_R1/2.fastq.gz: adapter trimmed fastq files.
* {prefix}_mkdup.bam: all alignments, with duplicates marked.
* {prefix}_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.
* {prefix}_se.bed: useful filtered alignments in BED format, Tn5 shifted.
* {prefix}_pe.bed: useful filtered alignments in BEDPE format,  Tn5 shifted, the 2nd and 3rd columns indicate the fragment start and end coordinates on genome. It will be used for macs2 peak calling.
* {prefix}.bw: bigwig file converted from {prefix}_filtered.bam, can be upload to genome browser for visualization.
* macs2: output of macs2, see [here](https://github.com/taoliu/MACS#output-files). Only broad peaks will be called by default. In addition, {prefix}_broad_filtered.bed is the peaks file with blacklist filtered.
* fastqc: the report(s) of fastqc
* logs: running logs

> No bam files and trimmed fastq files will be generated with chromap runs.



-----

## ChIPseq.sh

This script performs quality control on fastq files and aligns reads to a reference genome using either BWA or chromap (bowtie2 for CUT&RUN), depending on the species specified with -g or the index provided with -i. It then converts alignments to filtered BAM/BED and bigwig formats but does NOT perform peak calling.

> This script works for both ChIPseq and CUT&RUN.

#### Input
Paired-end fastq files with **_R1/2** suffix, i.e. test_R1.fastq.gz, test_R2.fastq.gz 
Or single-end fastq file with `-s`.

#### Options

help message can be shown by `ChIPseq.sh -h`

```shell
Usage: ChIPseq.sh <options> -g <hg38|hg19|mm10> <reads1>|..<reads2> 
  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -x [str] Custom BWA index PATH (valid only without -g option)
    -z [str] Custom chromosome size table (valid only without -g option)
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -s Single-end mod (Paired-end default)
    -n Nextera adapters (Truseq default)
    -a Use BWA aln algorithm (BWA mem default)
    -u CUR&RUN mode, will be paired-end mode and use bowtie2 aligner with --dovetail
    -b [str] Custom Bowtie2 index PATH  (valid only with -b option)
    -c Using chromap to process FASTQ instead of canonical bowtie2/bwa
    -i [str] Custom chromap genome index (valid only with -c option)
    -r [str] Custom chromap genome reference (valid only with -c option)
```

#### Example

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ChIPseq.sh
chmod 755 ChIPseq.sh
./ChIPseq.sh -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```
Alternatively, you may use chromap aligner to speed up the processing :

```shell
./ChIPseq.sh -c -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

####  Output

All results will be store in current (./) directory.

* {prefix}_trimmed_R1/2.fastq.gz: adapter trimmed fastq files.
* {prefix}_mkdup.bam: all alignments, with duplicates marked.
* {prefix}_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.
* {prefix}_se.bed: useful filtered alignments in BED format.
* {prefix}_pe.bed: useful filtered alignments in BEDPE format, the 2nd and 3rd columns indicate the fragment start and end coordinates on genome.
* {prefix}.bw: bigwig file converted from {prefix}_filtered.bam, can be upload to genome browser for visualization.
* fastqc: the report(s) of fastqc
* logs: running logs

> No bam files and trimmed fastq files will be generated with chromap runs.


#### Peak calling
> NOTE:
> this pipeline does NOT call peaks, you may run it manually.
> Input is highly recommended for peak calling, process input fastq files with this pipeline with same parameter(s).

test_pe.bed (and input_pe.bed) can be used for macs2 peak calling in BEDPE mode:

```shell
macs2 callpeak -t test_pe.bed -c input_pe.bed -f BEDPE -g hs -n test
```

> --broad is recommended for histone modifications.

test_se.bed and test_filtered.bam can also be used in BED or BAM mode of macs2.

See more about [MACS2](https://github.com/taoliu/MACS) (for TFs peak calling) and [SICER](https://zanglab.github.io/SICER2/) (for Histone Mods peak calling).



-----


## RNAseq.sh

This script performs quality control on fastq files, aligning reads to either the reference genome or transcriptome using STAR, based on the species selected via -s, or using the specified index and GTF files via -i and -g. Reads are counted using featureCounts, and differential expression analysis is conducted using DESeq2 to discover differentially expressed genes.

#### Input

* Paired-end fastq files with _R1/2.fastq.gz or _1/2.fq.gz suffix all together in a directory, i.e.

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

And a text file with meta information. i.e.

```shell
Sample  Group
cond1_rep1  group1
cond1_rep2  group1
cond2_rep1  group2
cond2_rep2  group2
....
```

You can use the script to scan fastq files and generate the text file:

```shell
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/RNAseq.sh
chmod 755 RNAseq.sh 
./RNAseq.sh -p /path/to/directory/contains/fastq/
```
Then the meta.txt will be created and opened with VIM. Sample column (1st) should have been filled, edit the text by adding the group information on the 2nd column.

> NOTE:
> Provide **the PATH of the DIRECTORY** which contains fastq to the scripts, DO NOT give the path of fastq files directly!

#### Options

help message can be shown by `RNAseq.sh -h`

```shell
Usage: RNAseq.sh <options> -m meta.txt </PATH/to/fastq/> 
  Options:
    -m [str] /PATH/to/meta.txt
    -s [str] species <hg|mm> hg=hg38, mm=mm10
    -i [str] Custom STAR index PATH
    -g [str] Custom Reference GTF transcripts PATH
    -t [int] Threads (1 default)
    -p prepare meta.txt
    -n Nextera adapters (Truseq default)
    -h Print this help message
```

#### Example

```shell
./RNAseq.sh -s hg -m meta.txt -t 24 /path/to/directory/contains/fastq/
```

Alternatively, you may use your custom genome and transcriptome reference:

```shell
./RNAseq.sh -i /path/to/STAR/index -g /path/to/GTF -m meta.txt -t 24 /path/to/directory/contains/fastq/
```

####  Output

All results will be store in current (./) directory.

* TRIMMED/{prefix}_R1/2_trimmed.gz: adapter trimmed fastq files.
* BAM/{prefix}.bam: STAR output, accepted alignments.
* BAM/SJ.out.tab: STAR output, splice junctions.
* featureCounts.txt: featureCounts output, raw fragments count.
* Allgene_TPM.csv: Transcripts-Per-Million values of all genes in the reference.
* DESeq2_A_vs_B.csv: Pairwise DESeq2 test results for the groups specified in *meta.txt*, with gene TPM values.
* logs: running logs and fastqc reports.

> NOTE:
Sample names in meta.txt have to match the featureCounts output exactly, check your text or use this script to create it.
This script automatically performs pairwise differential gene expression (DEG) analysis when there are more than two groups. A online tool might be useful: [iDEP](http://bioinformatics.sdstate.edu/idep/).



----

## adapt_trim.sh

This script is separated from ChIPseq.sh, it trims adapter sequences from fastq files with cutadapt@python3.

#### Input
* Paired-end fastq files or single-end with -s, i.e. test_R1.fastq.gz test_R2.fastq.gz

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

* {prefix}_R1/2_trimmed.gz: adapter trimmed fastq files.

  

------
## trans_assemble.sh

This script QC fastq files and aligns reads to hg19/GRCh37(depends on the index and GTF provided) using HISAT2. *De novo* transcripts assembly will be performed by stingtie.

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

This script is a re-written wrapper of [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR).  
> Requirements
stringtie, cuffdiff, plar, CPC2 and HMMER

#### Input
* BAM files of aligned RNAseq data, and coresponding GTFs.

#### Usage

```shell
./PLAR.sh hornet.bam <output dir> <prefix_sample1,prefix_sample2...> <strand:rf|fr|un> <sample1_rep1.gtf> <sample1_rep2.gtf> <sample2_rep1.gtf> ... <sample1_rep1.bam,sample1_rep2.bam,sample2_rep1..>
```

> NOTE:
Edit the script and modify $plar_path, $cpc2_path.
Additional annotation files required in $plar_path, see [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) for more details.

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

## CRISPRlib.sh

This script uses cutadapt trimming the input fastq files to get the sgRNA sequences (20nt) according to the adaptor sequence on the 5' end next to the sgRNA, and then aligns these sequences to bowtie index build with the reference sgRNA library.  The Addgene library 110066, 160129 and 162256 have been prebuilt and can be directly assigned by  -l.

#### Input

* The fastq file containing sgRNA sequences.

#### Usage

help message can be shown by `CRISPRlib.sh -h`

```shell
    Usage: CRISPRlib.sh <options> <reads_clean.fq.gz>

    ### INPUT: fastq files ###
    This script will trim the input fastq to 20nt after the given sequence with cutadapt, and align the trimmed reads to the reference library build with Bowtie1, depending on the library selection passed by -l or the index and adapter sequence passed by -i and -a,
    then statisticize each sequence's frequency, and all results will be store in current (./) directory.
    ### python3/cutadapt/bowtie1/samtools required ###

    Options:
    -l [str] library selection <110066|160129|162256>
    -i [str] Custom bowtie index PATH
    -a [str] Custom adapter sequence
    -p [str] Prefix of output
    -n [int] Threads (1 default)
    -h Print this help message
```

#### Example

```shell
CRISPRlib.sh -n 12 -l 110066 -p test lib_R2.fastq.gz
```

Alternatively, you may build your custom sgRNA library index and give the path to the script along with your adaptor sequence:

1. First make sure your sgRNA library is in [FASTA](https://www.ncbi.nlm.nih.gov/genbank/fastaformat/) format, such as:

```shell
>sgZC3H12A_5
GGGCAGCGACCTGAGACCAG
>sgZC3H12A_6
GGAGTGGAAGCGCTTCATCG
...
```

Here is an example to format a csv file (let's say sgRNA name is in 1st column and sequence in 2nd column) to FASTA:

```shell
awk -F',' '{print ">"$1"\n"$2}' sample.csv > sample.fa
```

2. Then build bowite1 index using this FASTA:

```shell
bowtie-build --threads 12 sample.fa sgLib
```

This will generate a few index files with the prefix "sgLib" 

3. Now you can run CRISPRlib.sh with your own library and adaptor sequence.

```shell
CRISPRlib.sh -n 12 -i /path/to/index/with/prefix -a YOURADAPTORSEQUENCE -p test lib_R2.fastq.gz
```

> The adaptor sequence is 5'-end next to the sgRNA, depends on the vector used. 8nt or more would be recommended.

> Note: edit line53 to 58 or insert new options in the script if you want to automatically assign your own index path and adaptor sequences with -l.

#### Output

All results will be store in current (./) directory.

* {prifix}.tr.fq.gz: trimmed fastq containing the sgRNA only (20bp)
* {prefix}.sam: bowtie1 alignment
* {prefix}.srt.bam:  bowtie1 alignment in sorted and indexed BAM format 
* {prefix}.srt.bam.bai:  bowtie1 alignment index
* {pre}.log: trimming and alignment summary
* {pre}.counts.tsv: a table delimited text with sgRNA name (1st column) and its count number (2nd column)
* {pre}.table.tsv: insert the actual sequence between the name and count compared to counts.tsv

------
Author [@zhaoshuoxp](https://github.com/zhaoshuoxp) 
April 24 2024
