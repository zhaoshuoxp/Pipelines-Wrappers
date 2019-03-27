# ATACseq, ChIPseq, RNAseq pipeline wrappers amd more
-----
This repository has the following combined shell/awk/python/R scripts which can be used for High-throughput sequecning data analysis.

> * ATACseq.sh: bulk ATACseq pipeline, from fastq to open chromatin regions.
> * ChIPseq.sh: ChIPseq pipeline, from fastq to peak calling step.
> * RNAseq.sh: bulk ATACseq pipeline, from fastq to differential expression genes.
> * adapt_trim.sh: adapter trimming function, seperated from the above pipelines.
> * cisVar.sh: pipeline wrapper of [cisVar](https://github.com/TheFraserLab/cisVar).
> * GATK_HF.sh: variants calling by [GATK](https://software.broadinstitute.org/gatk/), from fastq to vcf.
> * trans_assemble.sh: *de novo* transcript assembly, from fatsq to GTF.
> * PLAR.sh: *de novo* [PLAR](http://www.weizmann.ac.il/Biological_Regulation/IgorUlitsky/PLAR) lncRNA discovery pipeline wrapper.
> * rRNA_dep.sh: ribosomal RNA depletion from fastq files.

-----
