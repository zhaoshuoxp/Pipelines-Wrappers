#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
which cutadapt &>/dev/null || { echo "cutadapt not found, install by pip3!"; exit 1; }
which python3 &>/dev/null || { echo "python3 not found, install python3!"; exit 1; }
which bowtie2 &>/dev/null || { echo "bowtie2 not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which samtools &>/dev/null || { echo "samtools not found!"; exit 1; }
which bedtools &>/dev/null || { echo "bedtools not found!"; exit 1; }
which macs2 &>/dev/null || { echo "macs2 > 2.1.1 not found!"; exit 1; }
which bedGraphToBigWig &>/dev/null || { echo "bedGraphToBigWig not found!"; exit 1; }
which bedItemOverlapCount &>/dev/null || { echo "bedItemOverlapCount not found!"; exit 1; }

#### DEFAULT CONFIGURATION ###
# default Paired-end mod
mod='pe'
# default Nextera adapters
aA='CTGTCTCTTATACACATCT'
gG='AGATGTGTATAAGAGACAG'
# default 1 core to run
threads=1
# Bowtie2 index
bw2index_hg19='/home/quanyi/genome/hg19/Bowtie2index/hg19.bowtie2'

# help message
help(){
	cat <<-EOF
  Usage: ATAC.sh <options> <reads1>|..<reads2> 

  ### INPUT: Paired-end fastq files with _R1/2 extension ###
  This script QC fastq files and align reads to hg19/GRCh37 using Bowtie2, 
  convert to filtered BAM/BED and bigwig format,
  then call peaks with MACS2 in BEDPE mode after Tn5 shifting,
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/bowtie2/samtools/bedtools/bedGraphToBigWig/bedItemOverlapCount/macs2(>=2.1.1) required ###

  Options:
    -i [str] Bowtie2 index PATH
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -s Single-end mod (DO NOT recommend, Paired-end default)
    -h Print this help message
EOF
	exit 0
}

QC_mapping(){
	if [ $1 = 'se' ];then
		# single-end CMD
		# FastQC 
		fastqc -f fastq -t $threads -o fastqc $3 
		# Nextera adapter trimming
		cutadapt -m 30 -j $threads -a $aA -g $gG -o ${2}_trimmed.fastq.gz $3 > ./logs/${2}_cutadapt.log
		# Bowtie2 align
		bowtie2 -X 2000 --local --mm -p $threads -x $bw2index_hg19 -U ${2}_trimmed.fastq.gz -S ${2}.sam
		echo 'Bowtie2 mapping summary:' > ./logs/${2}_align.log
		tail -n 15 nohup.out >> ./logs/${2}_align.log
	else
		# paired-end CMD
		# FastQC
		fastqc -f fastq -t $threads -o fastqc $3 $4
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${2}_trimmed_R1.fastq.gz -p ${2}_trimmed_R2.fastq.gz $3 $4 > ./logs/${2}_cutadapt.log
		# Bowtie2 align
		bowtie2 -X 2000 --local --mm -p $threads -x $bw2index_hg19 -1 ${2}_trimmed_R1.fastq.gz -2 ${2}_trimmed_R2.fastq.gz -S ${2}.sam
		echo 'Bowtie2 mapping summary:' > ./logs/${2}_align.log
		tail -n 15 nohup.out >> ./logs/${2}_align.log
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
		# remove chrM alignments
		samtools index $threads ${1}_rm.bam 
		samtools idxstats ${1}_rm.bam | cut -f 1 |grep -v M | xargs samtools view -b -@ $threads -o ${1}_chrM.bam ${1}_rm.bam
		# filter out unmapped/failedQC/secondary/duplicates alignments
		samtools view -@ $threads -f 2 -F 1796 -b -o ${1}_filtered.bam ${1}_chrM.bam
		echo >> ./logs/${1}_align.log
		echo 'flagstat after filter:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_filtered.bam >> ./logs/${1}_align.log
		# clean
		rm ${1}_rm.bam ${1}_rm.bam.bai 
	# paired-end CMD
	else
		# download picard.jar for PE duplicates removal
		wget https://github.com/broadinstitute/picard/releases/download/2.18.27/picard.jar
		# mark duplicates
		java -jar picard.jar MarkDuplicates INPUT=${1}_srt.bam OUTPUT=${1}_mkdup.bam METRICS_FILE=./logs/${1}_dup.log REMOVE_DUPLICATES=false
		echo 'flagstat after mkdup:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_mkdup.bam >> ./logs/${1}_align.log
		# remove chrM alignments
		samtools index -@ $threads ${1}_mkdup.bam 
		samtools idxstats ${1}_mkdup.bam | cut -f 1 |grep -v M | xargs samtools view -b -@ $threads -o ${1}_chrM.bam ${1}_mkdup.bam
		# filter our unmapped/failedQC/unpaired/duplicates/secondary alignments
		samtools view -@ $threads -f 2 -F 1804 -b -o ${1}_filtered.bam ${1}_mkdup.bam
		echo >> ./logs/${1}_align.log
		echo 'flagstat after filter:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_filtered.bam >> ./logs/${1}_align.log
		# sort bam by query name for bedpe 
		samtools sort -n -@ $threads -o ${1}.bam2 ${1}_filtered.bam
		# bam2bedpe
		bamToBed -bedpe -i ${1}.bam2 > ${1}.bedpe
		# clean
		rm ${1}_srt.bam ${1}.bam2 picard.jar 
	fi
	rm ${1}.bam ${1}.sam ${1}_chrM.bam 
}

# Peak calling with MACS2 >v2.1.1
peak_calling(){
	if [ $1 = 'se' ];then
		# Tn5 shift in SE mode
		awk -F $'\t' 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4}else if($6=="-"){${3}=${3}-5} print $0}' ${2}_se.bed > ${2}_shift_se.bed
		# broad peak calling
		cd macs2
		macs2 callpeak -t ../${2}_shift_se.bed -g hs -n ${2} -f BED --keep-dup all --broad -B --SPMR --nomodel --shift -37 --extsize 73
	else
		# Tn5 shift in PE mode
		awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){if($2>=5){print $1,$2-5,$6-5}else if($6>5){print $1,0,$6-5}}}' ${2}.bedpe > ${2}_shift.bed
		# broad peak calling
		cd macs2
		macs2 callpeak -t ../${2}_shift.bed -g hs -n ${2} -f BEDPE --keep-dup all --broad -B --SPMR 
	fi
}

# Blacklist filter 
blacklist_filter(){
	curl -s http://hgdownload.cse.ucsc.edu/goldenPath/hg19/encodeDCC/wgEncodeMapability/wgEncodeDacMapabilityConsensusExcludable.bed.gz | gunzip -c > bklt_hg19
	intersectBed -v -a ${1}_peaks.broadPeak -b bklt_hg19 > ${1}_broad_filtered.bed
	rm bklt_hg19
	cd ..
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "st:hi:p:" arg
do
	case $arg in
		# Bowtie2 index PATH
		i) bw2index_hg19=$OPTARG;;
		t) threads=$OPTARG;;
		# single-end mod
		s) mod='se';;
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

	QC_mapping $mod $prefix $1 $2

	sam_bam_bed $prefix $mod

	bam2bigwig $prefix ${prefix}_filtered.bam
	
	if [ ! -d macs2 ];then 
		mkdir macs2
	fi 

	peak_calling $mod $prefix

	blacklist_filter $prefix
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