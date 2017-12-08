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
picard=/home/quanyi/app/picard.jar
mkdir logs
mkdir fastqc

# fastqc control
fastqc -f fastq -o fastqc $1 
fastqc -f fastq -o fastqc $2

# cutadapt to trim adaptors
cutadapt -f fastq -m 30 -a CTGTCTCTTATACACATCT -A CTGTCTCTTATACACATCT -g AGATGTGTATAAGAGACAG -G AGATGTGTATAAGAGACAG -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > ./logs/$3_cutadapt.log

# bowtie2 aligment #up to: 4 aligment/insert 2000bp
bowtie2 -k 4 -X 2000 --local --mm -p 4 -x $bowtie2index_hg19 -1 $3_R1_trimmed.gz -2 $3_R2_trimmed.gz -S $3.sam
echo 'Bowtie2 mapping summary:' > ./logs/$3_align.log
tail -n 15 nohup.out >> ./logs/$3_align.log

# samtools to BAM/sort
samtools view -b $3.sam -o $3.bam
samtools sort $3.bam -o $3_srt.bam

# picard mark duplaicates
java -jar $picard MarkDuplicates INPUT=$3_srt.bam OUTPUT=$3_mkdup.bam METRICS_FILE=./logs/$3_dup.log

echo >> ./logs/$3_align.log
echo 'flagstat after markdup:' >> ./logs/$3_align.log
samtools flagstat $3_mkdup.bam >> ./logs/$3_align.log
samtools index $3_mkdup.bam 

# remove chrM reads
samtools idxstats $3_mkdup.bam | cut -f 1 |grep -v M | xargs samtools view -b -o $3_chrM.bam $3_mkdup.bam
samtools index $3_chrM.bam

# remove unpaired/unmapped/duplicates/not primary alignments
samtools view -f 2 -F 1804 -b -o $3_filtered.bam $3_chrM.bam
echo >> ./logs/$3_align.log
echo 'flagstat after filter:' >> ./logs/$3_align.log
samtools flagstat $3_filtered.bam >> ./logs/$3_align.log

# clean
rm $3.sam $3.bam $3_srt.bam $3_chrM.bam $3_chrM.bai $3_R1_trimmed.gz $3_R2_trimmed.gz
