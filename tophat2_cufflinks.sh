#!/bin/bash
#####################################
# Usage:	$1=reads1/file/to/path	#
#			$2=reads2/file/to/path	#
#			$3=output_file_prefix	#
#####################################

# path to reads
READS1=$1
READS2=$2
NAME=$3
mkdir $3
mkdir $3/fastqc

# fastqc control
fastqc -f fastq -o $3/fastqc $1 
fastqc -f fastq -o $3/fastqc $2

# cutadapt--trim adaptors
cutadapt -f fastq -m 20 -a CTGTCTCTTATA -A CTGTCTCTTATA -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > ./$3/$3_cutadapt.log

# tophat--mapping
tophat2 -r 150 -p 8 --transcriptome-index=$transcriptome_bowtie2_dm6 -o $3 $bowtie2index_dm6 $3_R1_trimmed.gz $3_R2_trimmed.gz

# cufflinks--call transcripts
cufflinks -o $3 -p 8 -g $gtf_dm6 -L $3 -u --total-hits-norm ./$3/accepted_hits.bam
mv ./$3/transcripts.gtf ./$3/$3.gtf 

# clean
rm $3_R1_trimmed.gz $3_R2_trimmed.gz

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################