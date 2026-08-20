# Single cell and bulk sequencing pipeline wrappers

-----
This repository has the following combined shell/awk/python/R scripts which can be used for High-throughput sequecning data analysis.

* [cellranger.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#cellrangersh): Run CellRanger in batches to process 10x Genomics single-cell fastq files.
 * [ATACseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#atacseqsh): bulk ATACseq and CUT&TAG pipeline, from fastq to open chromatin regions/peaks.
 * [ChIPseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#chipseqsh): ChIPseq pipeline, from fastq to peak calling step.
 * [RNAseq.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rnaseqsh): bulk RNAseq pipeline, from fastq to differentially expressed genes.
 * [chrombpnet.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#chrombpnetsh): run chrombpnet workflow, from fragment files to nobias models.
 * [LDlookup.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#ldlookupsh): find Linkage Disequilibrium (LD) proxy SNPs for specific populations.
 * [adapt_trim.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#adapt_trimsh): adapter trimming function, seperated from the above pipelines.
 * [cisVar.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#cisvarsh): pipeline wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).
 * [TPS_assemble.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#tps_assemblesh): *de novo* transcript assembly, from fastq to GTF.
 * [PLAR.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#plarsh): *de novo* [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) lncRNA discovery pipeline wrapper.
 * [rRNA_dep.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#rrna_depsh): ribosomal RNA depletion from fastq files.
 * [CRISPRlib.sh](https://github.com/zhaoshuoxp/Pipelines-Wrappers#crisprlibsh): mapping CRISPR sgRNA library, from fastq to tables.

> Requirements:
> Python3, cutadapt, macs2, R, DESeq2, featureCounts, bowtie1/2, bwa,STAR, fastqc, samtools, bedtools, deeptools, chrombpnet conda env, cellranger(10 by default), cellranger-atac (2.2) and cellranger-arc

[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu) [![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

-----

## cellranger.sh

This script runs Cell Ranger count or aggr for RNA, ATAC, and multiome data.

For RNA/ATAC count, omitting `--id` and `--sample` enables automatic mode. The script recursively scans the FASTQ root, infers sample prefixes from standard filenames such as:

```
<prefix>_S*_L*_R?_*.fastq.gz
```

FASTQs may be stored together or in sample subdirectories. All lanes with the same exact prefix are passed to one count invocation without merging files. Samples are processed sequentially.

Multiome count accepts either `--libraries CSV` or separate `--gex_path` and `--atac_path` roots.

Aggr accepts an existing CSV with `--csv`. If omitted, complete count directories under `COUNT_ROOT` are discovered and an aggregation CSV is generated automatically. Normalization defaults to `none`.

### Usage

```
cellranger.sh -m <rna|atac> -g GENOME [options] FASTQ_ROOT
cellranger.sh -m <rna|atac> -g GENOME --sample PREFIX [--id ID] FASTQ_DIRS
cellranger.sh -m multiome -g GENOME --libraries CSV [options]
cellranger.sh -m multiome -g GENOME --gex_path DIR --atac_path DIR [options]
cellranger.sh -m <rna|atac|multiome> -a [--csv CSV | COUNT_ROOT] [options]
```

Use `cellranger.sh -h` for the complete option list.

Important options:

```
-m MODE              rna, atac, or multiome
-g GENOME            hg38, mm10, or mm39
-x PATH              Custom reference; overrides -g
-t N                 Cores; default: 20
-r GB                Memory; default: 200
-u VERSION           RNA Cell Ranger version: 7, 8, 9, or 10
--id ID              Optional run ID; generated when omitted
--sample PREFIX      Explicit RNA/ATAC FASTQ prefix
--libraries CSV      Multiome libraries CSV
--gex_path DIR       Multiome GEX FASTQ root
--atac_path DIR      Multiome ATAC FASTQ root
-a                    Run aggr
-c, --csv CSV        Aggregation CSV
--normalize MODE     none, mapped for RNA, or depth for ATAC/multiome
--output-dir DIR     Output path or automatic-mode output root
--dry                Validate and preview commands
```

### Examples

Automatic ATAC count:

```
./cellranger.sh -m atac -g mm10 -t 12 -r 96 /path/to/atac_fastqs
```

Explicit RNA count with multiple lanes:

```
./cellranger.sh -m rna -g mm10 \
  --sample MySample \
  /path/to/lane1,/path/to/lane2
```

Automatic multiome count:

```
./cellranger.sh -m multiome -g mm10 \
  --gex_path /path/to/GEX \
  --atac_path /path/to/ATAC
```

ATAC aggr without an existing CSV:

```
./cellranger.sh -m atac -g mm10 -a \
  --normalize none \
  /path/to/atac_count_outputs
```

### Output and restart behavior

If `--output-dir` is omitted, the current directory is used as the output root. Each result is stored under `RUN_ID/`.

Generated multiome libraries CSV and aggregation CSV files are written directly to the current directory.

Automatic status is appended to:

```
cellranger_auto_rna/status.log
cellranger_auto_atac/status.log
cellranger_auto_multiome/status.log
```

Monitor with:

```
tail -f cellranger_auto_atac/status.log
```

Re-running the same command skips complete outputs and resumes incomplete pipestances using deterministic run IDs.

> Sequencing-depth normalization is disabled by default. For large RNA atlases, combine count matrices in Seurat or Scanpy instead of using Cell Ranger RNA aggr.



-----

## ATACseq.sh

This script checks quality controls of fastq files, then aligns reads to the specified reference genome using Bowtie2 or chromap, depending on the selected species passed by -g or the provided index and other necessary files specified by -i, -b, and -c. It converts the alignments to filtered BAM/BED and bigwig formats, and subsequently identifies peaks using MACS2 in BED mode following Tn5 shifting.

#### Input
Paired-end fastq files with **_R1/2** suffix, i.e. test_R1.fastq.gz, test_R2.fastq.gz 
> Single-end sequencing data is also supported with -s, but it is not recommended.

#### Options

help message can be shown by `ATACseq.sh -h`

```
 Usage: ATAC.sh <options> <reads1>|<reads2>

  ### INPUT: Paired-end fastq files ###
  This script will QC fastq files and align reads to reference genome build with Bowtie2 or chromap, depending on the species passed by -g or the index and other required files passed by -i, -b and -c, convert alignments to filtered BAM/BED and bigwig, then call peaks with MACS2 in BED mode after Tn5 shifting.
  ### python3/cutadapt/fastqc/bowtie2/samtools/bedtools/deeptools/macs2 required ###

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

```
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ATACseq.sh
chmod 755 ATACseq.sh
./ATACseq.sh -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

Alternatively, you may use [chromap](https://github.com/haowenz/chromap) aligner to speed up the processing:

```
./ATACseq.sh -c -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```

> chromap is recomended for most cases, as it is ultra fast and outputs a little bit less (<5%) aligned reads.

####  Output

All results will be store in current (./) directory.

* {prefix}_trimmed_R1/2.fastq.gz: adapter trimmed fastq files.
* {prefix}_mkdup.bam: all alignments, with duplicates marked.
* {prefix}_filtered.bam: useful filtered alignments; duplicates, unpaired, unmapped, low-quality, secondary, chrM reads removed.
* {prefix}_se.bed: useful filtered alignments in BED format, Tn5 shifted. It will be used for macs2 peak calling.
* {prefix}_pe.bed: useful filtered alignments in BEDPE format,  Tn5 shifted, the 2nd and 3rd columns indicate the fragment start and end coordinates on genome. 
* {prefix}.bw: bigwig file converted from {prefix}_filtered.bam, can be upload to genome browser for visualization.
* macs2: output of macs2, see [here](https://github.com/taoliu/MACS#output-files). Only narrow peaks will be called by default. In addition, {prefix}_filtered.bed is the peaks file with blacklist filtered.
* fastqc: the report(s) of fastqc
* logs: running logs

> No bam files and trimmed fastq files will be generated with chromap runs.



-----

## ChIPseq.sh

This script performs quality control on fastq files and aligns reads to a reference genome using either BWA or chromap (bowtie2 for CUT&RUN/CUT&TAG), depending on the species specified with -g or the index provided with -i. It then converts alignments to filtered BAM/BED and bigwig formats but does NOT perform peak calling.

> This script works for both ChIPseq and CUT&RUN/CUT&TAG.

#### Input
Paired-end fastq files with **_R1/2** suffix, i.e. test_R1.fastq.gz, test_R2.fastq.gz 
Or single-end fastq file with `-s`.

#### Options

help message can be shown by `ChIPseq.sh -h`

```
  Usage: ChIPseq.sh <options> <reads1>|<reads2>

  ### INPUT: Single-end or Paired-end fastq files ###
  This script will QC fastq files and align reads to reference genome with BWA or chromap (bowtie2 for CUT&RUN/TAG), depending on the species passed by -g or the index passed by -i, convert alignments to filtered BAM/BED and bigwig but DOES NOT call peaks.
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/bwa/samtools/bedtools/deeptools required ###

  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -x [str] Custom BWA index PATH (valid only without -g option)
    -z [str] Custom chromosome size table (valid only without -g option)
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -s Single-end mod (Paired-end default)
    -n Manually force Nextera adapters (overrides defaults)
    -a Use BWA aln algorithm (BWA mem default)
    -R CUT&RUN mode: Paired-end, Bowtie2, TruSeq adapters
    -T CUT&Tag mode: Paired-end, Bowtie2, Nextera adapters
    -b [str] Custom Bowtie2 index PATH  (valid only with -u/-T option)
    -c Using chromap to process FASTQ instead of canonical bowtie2/bwa
    -i [str] Custom chromap genome index (valid only with -c option)
    -r [str] Custom chromap genome reference (valid only with -c option)
    -h Print this help message

```

#### Example

```
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/ChIPseq.sh
chmod 755 ChIPseq.sh
./ChIPseq.sh -g hg19 -p test -t 24 /path/to/test_R1.fastq.gz /path/to/test_R2.fastq.gz
```
Alternatively, you may use chromap aligner to speed up the processing :

```
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

```
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

```
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

```
Sample  Group
cond1_rep1  group1
cond1_rep2  group1
cond2_rep1  group2
cond2_rep2  group2
....
```

You can use the script to scan fastq files and generate the text file:

```
wget https://raw.githubusercontent.com/zhaoshuoxp/Pipelines-Wrappers/master/RNAseq.sh
chmod 755 RNAseq.sh 
./RNAseq.sh -p /path/to/directory/contains/fastq/
```
Then the meta.txt will be created and opened with VIM. Sample column (1st) should have been filled, edit the text by adding the group information on the 2nd column.

> NOTE:
> Provide **the PATH of the DIRECTORY** which contains fastq to the scripts, DO NOT give the path of fastq files directly!

#### Options

help message can be shown by `RNAseq.sh -h`

```
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

```
./RNAseq.sh -s hg -m meta.txt -t 24 /path/to/directory/contains/fastq/
```

Alternatively, you may use your custom genome and transcriptome reference:

```
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



-----

## chrombpnet.sh

This script automates the complete ChromBPNet workflow, from peak calling to final model training. It is designed to handle multiple clusters/cell types in parallel and features an advanced, dynamic peak downsampling module to balance single-cell data.

#### Input

Prepare a text file (e.g., `input_list.txt`) with two columns (tab or space separated):

```
Endothelial   /path/to/fragments/Endothelial.tsv.gz
SMC           /path/to/fragments/SMC.tsv.gz
Fibroblast    /path/to/fragments/Fibroblast.tsv.gz
```

- **Column 1:** Cluster Name (used for output directory and file naming).
- **Column 2:** Absolute path to the fragment file (`.tsv.gz`).

------

#### Standard Run (Auto-merge Bias)

If you don't have a merged fragment file for the bias model, the script will create one automatically by merging all fragments in your input list. **⚠️ Warning:** Ensure your temp directory partition has enough disk space (2-3x total size of fragments) for sorting.

```
conda activate chrombpnet
./chrombpnet.sh -i input_list.txt
```

### Run with Advanced Peak Downsampling 

When merging peaks from multiple single-cell clusters, the total peak count can easily exceed the optimal range for ChromBPNet training (typically 150k - 300k). Enable the advanced downsampling engine (`-s`) to automatically evaluate mathematical inflection points (Robust Scaled Composite Score) and apply proportional quotas to protect rare cell populations.

```
# Enable downsampling with the default target of 300,000 peaks
./chrombpnet.sh -i input_list.txt -s

# Enable downsampling and set a custom target peak count (e.g., 150,000)
./chrombpnet.sh -i input_list.txt -s -t 150000
```

#### Run with Existing Bias File

If you already have a merged bias fragment file, specify it with `-b` to save time.

```
./chrombpnet.sh -i input_list.txt -b /path/to/merged_fragments.tsv.gz
```

#### Custom Genome / Configuration

You can override default paths (hg38) using flags:

```
./chrombpnet.sh \
    -i input_list.txt \
    -g /path/to/mm10.fa \
    -c /path/to/mm10.chrom.sizes \
    -l /path/to/mm10.blacklist.bed
```

------

#### Arguments

| Flag     | Description                                                  | Default                                              |
| -------- | ------------------------------------------------------------ | ---------------------------------------------------- |
| **`-i`** | **[Required]** Input file list                               | N/A                                                  |
| **`-s`** | Enable advanced downsampling (Robust Scaled Knee + Proportional Fallback) | Disabled                                             |
| **`-t`** | Target peak count for downsampling (requires `-s`)           | 300000                                               |
| **`-b`** | Bias fragments file (`.tsv.gz`)                              | Auto-merge from input list                           |
| **`-g`** | Genome FASTA                                                 | `/home/quanyiz/genome/hg38/hg38.fa`                  |
| **`-c`** | Chromosome sizes                                             | `/home/quanyiz/genome/hg38/hg38.chrom.sizes`         |
| **`-l`** | Blacklist BED                                                | `/home/quanyiz/genome/hg38/hg38-blacklist.v2.bed`    |
| **`-f`** | Fold 0 JSON                                                  | `/home/quanyiz/genome/bias_models/folds/fold_0.json` |
| **`-d`** | Folds directory                                              | `/home/quanyiz/genome/bias_models/folds`             |

Export to Sheets

------

#### Output Structure

All outputs are generated in the `./results` folder or the directory where the script is executed:

```
./
├── macs2/          # MACS2 peaks, logs, and final union peak files
├── bias/           # Global bias model output (bias.h5)
├── models/         # Final ChromBPNet models (organized per cluster/fold)
└── tmp/            # Intermediate files (can be safely deleted after a successful run)
```



------

## LDlookup.sh

A robust, efficient, and highly automated Bash script designed to find Linkage Disequilibrium (LD) proxy SNPs for specific populations using local, [high-coverage 1000 Genomes Project (1kGP) VCF data](http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/) in hg38.

By leveraging local `bcftools` and `PLINK1.9`, this script completely bypasses the rate limits and network latency of online APIs (like LDlink), making it ideal for processing large batches of SNPs in seconds.

Ensure the following tools are installed and available in your Linux/Unix environment `$PATH`:

1. **`bcftools`**: For rapid VCF region extraction and sample subsetting.
2. **`plink` (v1.9 recommended)**: For blazing-fast LD calculations.
3. **`curl` & `python3`**: Used strictly for the Ensembl API coordinate fetch (only required if your input contains rsIDs).
4. **1kGP hg38 Reference Data** http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/:
   - Chromosome-level `.vcf.gz` files and their indices (`.tbi` or `.csi`).
   - The population pedigree file: `20130606_g1k_3202_samples_ped_population.txt`.

------

#### Basic Usage

The most basic run only requires your input list (`-i`) and the target population code (`-p`):

```
chmod +x run_proxy_ld.sh
./LDlookup.sh -i target_snps.txt -p AFR
```

#### Custom Parameters

```
./LDlookup.sh -i target_snps.txt -p EUR -m 0.05 -w 250000 -r 0.8 -o EUR_Proxies_R2_08.tsv
```

#### Command-Line Arguments

| Flag | Long Name   | Description                                                  | Default Value                       |
| ---- | ----------- | ------------------------------------------------------------ | ----------------------------------- |
| `-i` | `--input`   | **[Required]** Path to the input plain text file containing your SNPs. | None                                |
| `-p` | `--pop`     | **[Required]** Target population code (e.g., `AFR`, `EUR`, `GBR`). | None                                |
| `-v` | `--vcf-dir` | Directory path containing the 1kGP VCF files.                | `/path/to/1000Genome_hg38`          |
| `-d` | `--ped`     | Full path to the 1kGP PED sample population file.            | `.../20130606_g1k...txt`            |
| `-m` | `--maf`     | Minimum Allele Frequency threshold to filter out ultra-rare variants. | `0.001`                             |
| `-w` | `--window`  | LD search window size upstream/downstream (in base pairs).   | `500000` (+/- 500kb = 1Mb total)    |
| `-r` | `--r2`      | Minimum *R*2 threshold for the output proxies.               | `0` (Outputs all pairs in window)   |
| `-o` | `--out`     | Output TSV file name.                                        | `All_Local_Proxy_SNPs_Combined.tsv` |

------

#### Input Format

The input file (`-i`) should be a plain text file with one SNP per line. The script is highly fault-tolerant and accepts a mix of formats:

```
# Supported formats:
rs1234567
10:44256882:A:C
chr13:110209271:G:A
16:72182890
chr19:44911142
```

*Note: Even if you provide specific alleles (like `A:C`), the script will verify the actual alleles in the VCF to ensure PLINK runs successfully. However, **your exact original input string will be perfectly preserved** in the output table, ensuring seamless downstream merging with your original datasets.*

------

#### Output Format

Upon completion, a tab-separated values (TSV) file is generated (default: `All_Local_Proxy_SNPs_Combined.tsv`), which can be directly imported into R, Python, or Excel.

It contains the following columns:

1. **`CHR_A`**: Chromosome of the Lead SNP.
2. **`BP_A`**: Physical position of the Lead SNP (hg38).
3. **`TARGET_SNP`**: **Your exact original input ID** (crucial for `merge()`/`join` operations).
4. **`CHR_B`**: Chromosome of the Proxy SNP.
5. **`BP_B`**: Physical position of the Proxy SNP (hg38).
6. **`PROXY_SNP`**: The Proxy SNP ID found (formatted as `chr:pos:ref:alt`).
7. **`R2`**: Linkage Disequilibrium *R*2 metric.
8. **`POPULATION`**: The population code used for this calculation.



-----

## adapt_trim.sh

This script is separated from ChIPseq.sh, it trims adapter sequences from fastq files with cutadapt@python3.

#### Input
* Paired-end fastq files or single-end with -s, i.e. test_R1.fastq.gz test_R2.fastq.gz

#### Options

help message can be shown by `adapt_trim.sh -h`

```
Usage: adapt_trim.sh <options> <reads1>|..<reads2
    Options:
    	-p Prefix of output
      -t Threads (1 default)
      -s Single-end mod (Paired-end default)
      -n Nextera adapters (Truseq default)
      -h Print this help message
```

#### Example

```
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
## TPS_assemble.sh

This script QC fastq files and aligns reads to hg19/GRCh37(depends on the index and GTF provided) using HISAT2. *De novo* transcripts assembly will be performed by stingtie.

#### Input

* Paired-end fastq files.

#### Usage

```shell
./TPS_assemble.sh <reads1> <reads2> <prefix of output> <starnd: fr|rf|un>
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

This script automates the preprocessing, alignment, and counting of CRISPR sgRNA reads. It takes a raw sgRNA library list (TSV/CSV) to automatically build a Bowtie1 index on the fly. Then, it uses `cutadapt` to trim the input FASTQ files, extracting the 20nt sgRNA sequences based on a user-provided 5' adapter. Finally, it aligns these trimmed reads to the generated index and summarizes the read counts for each sgRNA.

### Input

1. **Reads:** A raw or clean FASTQ file (`.fq` or `.fq.gz`) containing the sequencing reads.
2. **Library:** A comma-separated (CSV) or tab-separated (TSV) file containing your sgRNA library details. The first column must be the sgRNA ID/Name, and the second column must be the sgRNA sequence.

#### Usage

You can display the help message by running `CRISPRlib.sh -h`.

```
Usage: CRISPRlib.sh <options> <reads_clean.fq.gz>

Options:
  -s [file] sgRNA sequence list (TSV or CSV format: ID <tab/comma> sequence) [Required]
  -a [str]  Custom adapter sequence [Required]
  -p [str]  Prefix of output
  -n [int]  Threads (12 default)
  -h        Print this help message
```

#### Example

First, ensure your sgRNA library is in a simple CSV or TSV format. For example (`sample.csv`):

```
sgZC3H12A_5,GGGCAGCGACCTGAGACCAG
sgZC3H12A_6,GGAGTGGAAGCGCTTCATCG
```

*(Note: You no longer need to manually convert this to FASTA or run `bowtie-build` yourself. The script handles this automatically.)*

Now you can run `CRISPRlib.sh` by simply providing your library file, adapter sequence, and FASTQ file:

Bash

```
CRISPRlib.sh -n 12 -s sample.csv -a YOURADAPTORSEQUENCE -p test_ lib_R2.fastq.gz
```

**Note on Adapters:** The adapter sequence (`-a`) should be the sequence on the 5'-end immediately preceding the sgRNA, which depends on the vector used. Providing 8nt or more is highly recommended for accurate trimming.

#### Output

All results will be stored in the current (`./`) directory.

- **`{prefix}ref.fa`**: The intermediate FASTA file automatically generated from your TSV/CSV.
- **`{prefix}ref_index.\*`**: The Bowtie1 index files built from the FASTA.
- **`{prefix}tr.fq.gz`**: Trimmed FASTQ containing only the extracted sgRNA sequences (20bp).
- **`{prefix}sam`**: Raw Bowtie1 alignment output.
- **`{prefix}bam`**: Bowtie1 alignment in BAM format.
- **`{prefix}srt.bam`**: Bowtie1 alignment in sorted and indexed BAM format.
- **`{prefix}srt.bam.bai`**: Bowtie1 alignment index.
- **`{prefix}log`**: Trimming and alignment summary and logs.
- **`{prefix}counts.tsv`**: A tab-delimited text table with the sgRNA name (1st column) and its total read count (2nd column).
- **`{prefix}table.tsv`**: A tab-delimited text table containing the read count (1st column), actual aligned sequence (2nd column), and sgRNA name (3rd column).

------
Author [@zhaoshuoxp](https://github.com/zhaoshuoxp) 
April 24 2024
