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
cutadapt -f fastq -m 30 -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > ./logs/$3_cutadapt.log

# bwa mem alignment
bwa mem -M -t 16 $bwaindex_hg19 $3_R1_trimmed.gz $3_R2_trimmed.gz > $3.sam 

# samtools convert/sort
samtools view -b -o $3.bam $3.sam 
samtools sort -o $3_srt.bam $3.bam
echo 'flagstat after alignment:' > ./logs/$3_align.log
samtools flagstat $3_srt.bam >> ./logs/$3_align.log
echo >> ./logs/$3_align.log

# picard markduplicates
java -jar $picard MarkDuplicates INPUT=$3_srt.bam OUTPUT=$3_rmdup.bam METRICS_FILE=./logs/$3_dup.log REMOVE_DUPLICATES=true
echo 'flagstat after rmdup:' >> ./logs/$3_align.log
samtools flagstat $3_rmdup.bam >> ./logs/$3_align.log
samtools index $3_rmdup.bam

# remove unpaired/unmapped/duplicates/not primary alignments
#samtools view -f 2 -F 1804 -b -o $3_filtered.bam $3_rmdup.bam
#echo >> ./logs/$3_align.log
#echo 'flagstat after filter:' >> ./logs/$3_align.log
#samtools flagstat $3_filtered.bam >> ./logs/$3_align.log

# convert to BED
# bamToBed -i $3_rmdup.bam > $3.bed

# BAM to BW 
bam2bigwig.sh $3_rmdup.bam

# peak calling macs2.0.9
# macs2 -t ChIP.bam -c input.bam -f BAMPE -g hs -n NAME --outdir macs2 -B 

rm $3.bam $3.sam $3_R1_trimmed.gz $3_R2_trimmed.gz

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################