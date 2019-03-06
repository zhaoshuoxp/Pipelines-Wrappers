#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
which cutadapt &>/dev/null || { echo "cutadapt not found, install by pip3!"; exit 1; }
which python3 &>/dev/null || { echo "python3 not found, install python3!"; exit 1; }
which bwa &>/dev/null || { echo "bwa not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which samtools &>/dev/null || { echo "samtools not found!"; exit 1; }
which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found!"; exit 1; }
which bedItemOverlapCount &>/dev/null || { echo "bedItemOverlapCount not found!"; exit 1; }

#### DEFAULT CONFIGURATION ###
# default Paired-end mod
mod='pe'
# default BWA mem algorithm
alg='mem'
# default TrueSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1
# BWA index
bwaindex_hg19='/home/quanyi/genome/hg19/BWAindex/hg19bwa'

# help message
help(){
	cat <<-EOF
	Usage: ChIPseq.sh <options> <reads1>|..<reads2> 
	QC fastq files and align to hg19/GRCh37 using BWA, convert to filtered BAM/BED and bigwig format. 
	-i BWA index PATH
	-p Prefix of output
	-t Threads (1 default)
	-s Single-end mod (Paired-end default)
	-a Use BWA aln algorithm (BWA mem default)
	-h Print this help message
EOF
	exit 0
}

### main pipeline ###

QC_mapping(){
	if [ $1 = 'se' ];then
		# single-end CMD
		# FastQC 
		fastqc -f fastq -t $threads -o fastqc $4 
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -g $gG -o ${3}_trimmed.fastq.gz $4 > ./logs/${3}_cutadapt.log
		# BWA aln 
		if [ $2 = 'aln' ];then
			bwa aln -t $threads -k 2 -l 18 $bwaindex_hg19 ${3}_trimmed.fastq.gz > ${3}.sai
			bwa samse $bwaindex_hg19 ${3}.sai ${3}_trimmed.fastq.gz  > ${3}.sam
			rm ${3}.sai
		# BWA mem
		else
			bwa mem -M -t $threads $bwaindex_hg19 ${3}_trimmed.fastq.gz > ${3}.sam
		fi
	else
		# paired-end CMD
		# FastQC
		fastqc -f fastq -t $threads -o fastqc $4 $5
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${3}_trimmed_R1.fastq.gz -p ${3}_trimmed_R2.fastq.gz $4 $5 > ./logs/${3}_cutadapt.log
		# BWA aln 
		if [ $2 = 'aln' ];then
			bwa aln -t $threads -k 2 -l 18 $bwaindex_hg19 ${3}_trimmed_R1.fastq.gz > ${3}_R1.sai
			bwa aln -t $threads -k 2 -l 18 $bwaindex_hg19 ${3}_trimmed_R2.fastq.gz > ${3}_R2.sai
			bwa sampe $bwaindex_hg19 ${3}_R1.sai ${3}_R2.sai ${3}_trimmed_R1.fastq.gz ${3}_trimmed_R2.fastq.gz  > ${3}.sam
			rm ${3}_R1.sai ${3}_R2.sai
		# BWA mem
		else
			bwa mem -M -t $threads $bwaindex_hg19 ${3}_trimmed_R1.fastq.gz ${3}_trimmed_R2.fastq.gz > ${3}.sam
		fi
	fi
}

# Convert BAM to BigWig
bam2bigwig(){
	# create bed from bam, requires bedtools bamToBed
	bamToBed -i $2 -split > ${1}_se.bed
	# get chromasome length
	curl -s ftp://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz | gunzip -c > hg19_len
	# create plus and minus strand bedgraph
	cat ${1}_se.bed | sort -k1,1 | bedItemOverlapCount hg19 -chromSize=hg19_len stdin | sort -k1,1 -k2,2n > accepted_hits.bedGraph
	# convert to bigwig
	bedGraphToBigWig accepted_hits.bedGraph hg19_len ${1}.bw
	# removing intermediery files
	rm accepted_hits.bedGraph hg19_len
}

# SAM2BAM and filtering to BED
sam2bed(){
	# sam2bam+sort
	samtools view -b -@ $threads -o ${1}.bam ${1}.sam 
	samtools sort -@ $threads -o ${1}_srt.bam ${1}.bam
	# single-end CMD
	if [ $2 = 'se' ];then
		# samtools rmdup module for SE duplicates removal
		samtools rmdup -s ${1}_srt.bam ${1}_rm.bam
		echo 'flagstat after rmdup:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_rm.bam >> ./logs/${1}_align.log
		# filter out unmapped/failedQC/secondary/duplicates alignments
		samtools view -@ $threads -f 2 -F 1796 -b -o ${1}_filtered.bam ${1}_rm.bam
		echo >> ./logs/${1}_align.log
		echo 'flagstat after filter:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_filtered.bam >> ./logs/${1}_align.log
		#clean
		rm ${1}.bam ${1}_rm.bam
	# paired-end CMD
	else
		# download picard.jar for PE duplicates removal
		wget https://github.com/broadinstitute/picard/releases/download/2.18.27/picard.jar
		# mark duplicates
		java -jar picard.jar MarkDuplicates INPUT=${1}_srt.bam OUTPUT=${1}_mkdup.bam METRICS_FILE=./logs/${1}_dup.log REMOVE_DUPLICATES=false
		echo 'flagstat after mkdup:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_mkdup.bam >> ./logs/${1}_align.log
		# filter our unmapped/failedQC/unpaired/duplicates/secondary alignments
		samtools view -@ $threads -f 2 -F 1804 -b -o ${1}_filtered.bam ${1}_mkdup.bam
		echo >> ./logs/${1}_align.log
		echo 'flagstat after filter:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_filtered.bam >> ./logs/${1}_align.log
		# sort bam by query name for bedpe 
		samtools sort -n -@ $threads -o ${1}.bam2 ${1}_filtered.bam
		# bam2bedpe
		bamToBed -bedpe -i ${1}.bam2 > ${1}.bedpe
		# bedpe to standard PE bed for macs2 peak calling (-f BEDPE)
		cut -f1,2,6 ${1}.bedpe > $3_pe.bed
		# clean
		rm ${1}.bam2 ${1}.bedpe
	fi
}

# no ARGs error
if [ $# -lt 1 ]
then
	help
	exit 1
fi

while getopts "sat:hi:p:" arg
do
	case $arg in
		# BWA index PATH
		i) bwaindex_hg19=$OPTARG;;
		t) threads=$OPTARG;;
		# single-end mod
		s) mod='se';;
		# BWA algorithm
		a) alg='aln';;
		p) prefix=$OPTARG;;
		h) help ;;
		?) help
			exit 1;;
	esac
done

# shift ARGs to reads
shift $(($OPTIND - 1))
# get prefix of output
if [ -z $prefix ];then
	echo "No -p <prefix> given, use file name as prefix"
	if [ $mod = 'se' ];then
		prefix=${1%.*}
	else
		prefix=${1%_R1*}
	fi
fi

# main
if [ ! -d logs ]
then 
mkdir logs
fi

if [ ! -d fastqc ]
then 
mkdir fastqc
fi 

QC_mapping $mod $alg $prefix $1 $2

sam2bed $prefix $mod

bam2bigwig $prefix ${prefix}_filtered.bam

# check running status
if [ $? -ne 0 ]; then
	help
	exit 1
else
	echo "Run succeed"
fi

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################