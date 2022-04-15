#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
requires=("cutadapt" "python3" "bwa" "fastqc" "samtools" "bedtools" "bamCoverage")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

#### DEFAULT CONFIGURATION ###
# default Paired-end mod
mod='pe'
# default BWA mem algorithm
alg='mem'
# default TruSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1
# BWA index
bwaindex_hg19='/mnt/date3/Project/zhaoqy/genome/hg19/BWAindex/hg19.fa'
bwaindex_hg38='/mnt/date3/Project/zhaoqy/genome/hg38/BWAindex/GRCh38.p13.genome.fa'
bwaindex_mm10='/mnt/date3/Project/zhaoqy/genome/mm10/BWAindex/mm10.fa'
# Picard
picard_url='https://github.com/broadinstitute/picard/releases/download/2.27.1/picard.jar'

# help message
help(){
	cat <<-EOF
  Usage: ChIPseq.sh <options> <reads1>|<reads2> 

  ### INPUT: Single-end or Paired-end fastq files with _R1/2 extension ###
  This script will QC fastq files and align reads to reference genome with BWA, depending on the species selection passed by -g or the index passed by -i, 
  convert to filtered BAM/BED and bigwig format but DOES NOT call peaks.
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/bwa/samtools/bedtools/deeptools required ###

  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -i [str] Custom BWA index PATH
    -p [str] Prefix of output
    -t [int] Threads (1 default)
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
			bwa aln -t $threads -k 2 -l 18 $bwaindex ${3}_trimmed.fastq.gz > ${3}.sai
			bwa samse $bwaindex ${3}.sai ${3}_trimmed.fastq.gz  > ${3}.sam
			rm ${3}.sai
		# BWA mem
		else
			bwa mem -M -t $threads $bwaindex ${3}_trimmed.fastq.gz > ${3}.sam
		fi
	else
		# paired-end CMD
		# FastQC
		fastqc -f fastq -t $threads -o fastqc $4 $5
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${3}_trimmed_R1.fastq.gz -p ${3}_trimmed_R2.fastq.gz $4 $5 > ./logs/${3}_cutadapt.log
		# BWA aln 
		if [ $2 = 'aln' ];then
			bwa aln -t $threads -k 2 -l 18 $bwaindex ${3}_trimmed_R1.fastq.gz > ${3}_R1.sai
			bwa aln -t $threads -k 2 -l 18 $bwaindex ${3}_trimmed_R2.fastq.gz > ${3}_R2.sai
			bwa sampe $bwaindex ${3}_R1.sai ${3}_R2.sai ${3}_trimmed_R1.fastq.gz ${3}_trimmed_R2.fastq.gz  > ${3}.sam
			rm ${3}_R1.sai ${3}_R2.sai
		# BWA mem
		else
			bwa mem -M -t $threads $bwaindex ${3}_trimmed_R1.fastq.gz ${3}_trimmed_R2.fastq.gz > ${3}.sam
		fi
	fi
}

# SAM2BAM and filtering to BED
sam_bam_bed(){
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
		# clean
		rm ${1}_rm.bam
	# paired-end CMD
	else
		# download picard.jar for PE duplicates removal
		wget $picard_url
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
		cut -f1,2,6 ${1}.bedpe > ${1}_pe.bed
		# clean
		rm ${1}_srt.bam ${1}.bam2 ${1}.bedpe picard.jar
	fi
	rm ${1}.bam ${1}.sam
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "sag:t:hi:p:" arg
do
	case $arg in
		g) if [ $OPTARG = "hg19" ]; then
			bwaindex=$bwaindex_hg19
		   elif [ $OPTARG = "hg38" ]; then
			bwaindex=$bwaindex_hg38
		   elif [ $OPTARG = "mm10" ]; then
			bwaindex=$bwaindex_mm10
		   else
			echo "Only support hg38, hg19 or mm10, or pass your custom genome index"
			exit 1
		   fi;;
		# BWA index PATH
		i) bwaindex=$OPTARG;;
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
main(){
	if [ ! -d logs ];then 
		mkdir logs
	fi

	if [ ! -d fastqc ];then 
		mkdir fastqc
	fi 
	
	QC_mapping $mod $alg $prefix $1 $2

	sam_bam_bed $prefix $mod

	# convert filtered BAM to CPM normalized bigWig with deeptools
	bamCoverage --binSize 10 -p $threads --normalizeUsing CPM -b ${prefix}_filtered.bam -o ${prefix}.bw
}

main $1 $2

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