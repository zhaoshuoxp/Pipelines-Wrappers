#!/bin/bash
#####################################
# Usage:	$1=reads1/file/to/path	#
#			$2=reads2/file/to/path	#
#			$3=output_file_name		#
#####################################

# path to reads
READS1=$1
READS2=$2
NAME=$3
#picard=/home/quanyi/app/picard.jar

# fastqc control
fastqc -f fastq $1 
fastqc -f fastq $2

# cutadapt to trim adaptors
cutadapt -f fastq -m 30 -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATT -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > $3_cutadapt.log

# bwa mem alignment
bwa mem -M -t 16 $bwaindex_hg19 $3_R1_trimmed.gz $3_R2_trimmed.gz > $3.sam 

# samtools convert/sort/rmdup
samtools view -@ 16 -Sb $3.sam -o $3.bam
samtools flagstat -@ 16 $3.bam > alignment.log
samtools sort -@ 16 $3.bam -o $3_srt.bam
echo "\n" >> alignment.log
samtools rmdup $3_srt.bam $3_rm.bam >> alignment.log

#java -jar $picard SortSam INPUT=$3.sam OUTPUT=$3_srt.bam SORT_ORDER=coordinate
#samtools flagstat -@ 16 $3.bam > alignment.log
#java -jar $picard MarkDuplicates INPUT=$3_srt.bam OUTPUT=$3_rm.bam METRICS_FILE=marked_dup_metrics.txt REMOVE_DUPLICATES=true  

# convert to BED
bamToBed -i $3_rm.bam > $3.bed

# peak calling macs2.0.9
# macs2 -t ChIP.bam -c input.bam -f bam -g hs -n $3 


################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################