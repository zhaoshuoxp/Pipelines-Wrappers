#!/bin/bash
#####################################
# Usage:	$1=reads1/file/to/path	#
#			$2=reads2/file/to/path	#
#			$3=output_file_prefix	#
#####################################
# check programs
which cutadapt &>/dev/null || { echo "cutadapt not found!"; exit 1; }
which bowtie2 &>/dev/null || { echo "bowtie2 not found!"; exit 1; }
which tophat2 &>/dev/null || { echo "tophat not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which cufflinks &>/dev/null || { echo "cufflinks not found!"; exit 1; }
#define sebnif path
sebnif=/home/quanyi/app/sebnif-1.3/sebnif
# path to reads
READS1=$1
READS2=$2
NAME=$3
mkdir -p $3/fastqc

# fastqc control
fastqc -f fastq -o $3/fastqc $1 
fastqc -f fastq -o $3/fastqc $2

# cutadapt--trim adaptors
cutadapt -f fastq -m 20 -a CTGTCTCTTATA -A CTGTCTCTTATA -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > ./$3/$3_cutadapt.log

# tophat--mapping
tophat2 -r 150 -p 8 --transcriptome-index=$transcriptome_bowtie2_hg19 -o $3 $bowtie2index_hg19 $3_R1_trimmed.gz $3_R2_trimmed.gz

# cufflinks for known genes
cufflinks -o $3 -p 8 -L $3 -u --total-hits-norm -g $gtf_hg19./$3/accepted_hits.bam
# cufflinks for de novo lncRNA (sebnif)
#cufflinks -o $3 -p 8 -L $3 -u --total-hits-norm ./$3/accepted_hits.bam
mv ./$3/transcripts.gtf ./$3/$3.gtf 

# convert GTF format to Sebnif compatable
#awk '{if (length($1<=5)){print $0}}' $3/$3.gtf > $3/$3_sebnif.gtf

# sebnif run
#$sebnif -k -X gmm -E -0.1,0.95 -o $3/$3_sebnif $3/$3_sebnif.gtf
#cp $3/$3_sebnif/novel.final.gtf $3/$3_lncRNA.gtf

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################