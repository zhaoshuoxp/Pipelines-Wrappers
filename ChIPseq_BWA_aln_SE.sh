#!/bin/bash
#####################################
# Usage:	$1=reads1/file/to/path	#
#			$2=output_file_prefix	#
#####################################
# check programs: 
which cutadapt &>/dev/null || { echo "cutadapt not found!"; exit 1; }
which bwa &>/dev/null || { echo "bwa not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which samtools &>/dev/null || { echo "samtools not found!"; exit 1; }
which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found!"; exit 1; }

# path to reads
READS1=$1
NAME=$2
#picard=/home/quanyi/app/picard.jar
threads=16

if [ ! -d logs ]
then 
mkdir logs
fi

if [ ! -d fastqc ]
then 
mkdir fastqc
fi 

# fastqc control
fastqc -f fastq -t $threads -o fastqc $1 

# cutadapt to trim adaptors
# python3 version required for -j
cutadapt -f fastq -m 20 -j $threads -a AGATCGGAAGAGC -g GCTCTTCCGATCT -o $2_trimmed.gz $1 > ./logs/$2_cutadapt.log

# bwa aln alignment
bwa aln -t $threads -k 2 -l 18 $bwaindex_hg19 $2_trimmed.gz > $2.sai
bwa samse $bwaindex_hg19 $2.sai $2_trimmed.gz > $2.sam

# samtools convert/sort
samtools view -b -@ $threads -o $2.bam $2.sam 
samtools sort -@ $threads -o $2_srt.bam $2.bam

# picard markduplicates
#java -jar $picard MarkDuplicates INPUT=$2_srt.bam OUTPUT=$2_rmdup.bam METRICS_FILE=./logs/$2_dup.log REMOVE_DUPLICATES=false
samtools rmdup -s $2_srt.bam $2_rmdup.bam

echo 'flagstat after mkdup:' >> ./logs/$2_align.log
samtools flagstat -@ $threads $2_rmdup.bam >> ./logs/$2_align.log

# remove unpaired/unmapped/duplicates/failedQC
# Do NOT do filtering when Pooled ChIPseq!!!
samtools index -@ $threads $2_rmdup.bam
samtools view -@ $threads -f 2 -F 1796 -b -o $2_filtered.bam $2_rmdup.bam

echo >> ./logs/$2_align.log
echo 'flagstat after filter:' >> ./logs/$2_align.log
samtools flagstat -@ $threads $2_filtered.bam >> ./logs/$2_align.log
bamToBed -i $2_filtered.bam > $2_filtered.bed

# BAM to BW 
#create SE bed from bam
bamToBed -i $2_rmdup.bam -split > $2.bed

#create plus and minus strand bedgraph
cat $2.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=$len_hg19 stdin | sort -k1,1 -k2,2n > $2.bedGraph

bedGraphToBigWig $2.bedGraph $len_hg19 $2.bw

# clean
rm $2_trimmed.gz $2.sai $2.sam $2.bam $2_srt.bam $2.bedGraph $2.bed

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################