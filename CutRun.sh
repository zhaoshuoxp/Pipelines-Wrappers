#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
requires=("cutadapt" "python3" "bowtie2" "fastqc" "samtools" "bedtools" "bamCoverage" "chromap")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

#### DEFAULT CONFIGURATION ###
# default Paired-end mod
aln='bowtie2'
# default BWA mem algorithm
# default TruSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1
# BWA index
bw2index_hg19='/nfs/baldar/quanyiz/genome/hg19/Bowtie2index/hg19.bowtie2'
bw2index_hg38='/nfs/baldar/quanyiz/genome/hg38/Bowtie2index/hg38bowtie2'
bw2index_mm10='/nfs/baldar/quanyiz/genome/mm10/Bowtie2index/mm10.bowtie2'
chromapindex_hg19='/nfs/baldar/quanyiz/genome/hg19/Chromapindex/hg19'
chromapindex_hg38='/nfs/baldar/quanyiz/genome/hg38/Chromapindex/hg38'
chromapindex_mm10='/nfs/baldar/quanyiz/genome/mm10/Chromapindex/mm10'
chromapref_hg19='/nfs/baldar/quanyiz/genome/hg19/GRCh37.p13.genome.fa'
chromapref_hg38='/nfs/baldar/quanyiz/genome/hg38/GRCh38.p14.genome.fa'
chromapref_mm10='/nfs/baldar/quanyiz/genome/mm10/GRCm38.p6.genome.fa'
chromszize_hg19='/nfs/baldar/quanyiz/genome/hg19/hg19.chrom.sizes'
chromszize_hg38='/nfs/baldar/quanyiz/genome/hg38/hg38.chrom.sizes'
chromszize_mm10='/nfs/baldar/quanyiz/genome/mm10/mm10.chrom.sizes'
# Picard
picard_url='https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar'

# help message
help(){
	cat <<-EOF
  Usage: CutRun.sh <options> <reads1>|<reads2> 

  ### INPUT: Paired-end fastq files ###
  This script will QC fastq files and align reads to reference genome with Bowtie2 or chromap, depending on the species passed by -g or the index passed by -i, 
  convert alignments to filtered BAM/BED and bigwig but DOES NOT call peaks.
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/bwa/samtools/bedtools/deeptools/chromap required ###

  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -x [str] Custom Bowtie2 index PATH
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -c Using chromap to process FASTQ instead of canonical bowtie2
    -i [str] Custom chromap genome index (only valid with -c option)
    -r [str] Custom chromap genome reference (only valid with -c option)
    -z [str] Custom chromosome size table
    -n Nextera adapters (Truseq default)
    -h Print this help message

EOF
	exit 0
}

### main pipeline ###

QC_mapping(){
	if [ $1 = 'chromap' ];then
		chromap -r $chromapref -x $chromapindex --trim-adapters --remove-pcr-duplicates -t $threads -1 $3 -2 $4 -o ${2}_pe.bed
		echo 'chromap mapping summary:' > ./logs/${2}_align.log
		tail -n 14 nohup.out >> ./logs/${2}_align.log
	else
		# FastQC
		fastqc -f fastq -t $threads -o fastqc $3 $4
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -A $aA -o ${2}_trimmed_R1.fastq.gz -p ${2}_trimmed_R2.fastq.gz $3 $4 > ./logs/${2}_cutadapt.log
		# Bowtie2 aln --dovetail for CUT&RUN
		bowtie2 --dovetail --threads $threads -X 1000 -x $bw2index -1 ${2}_trimmed_R1.fastq.gz -2 ${2}_trimmed_R2.fastq.gz -S ${2}.sam
		echo 'Bowtie2 mapping summary:' > ./logs/${2}_align.log
		tail -n 15 nohup.out >> ./logs/${2}_align.log
	fi
}

# SAM2BAM and filtering to BED
sam_bam_bed(){
	if [ $2 = 'chromap' ];then
		awk 'substr($1,1,3)=="chr"' ${1}_pe.bed > ${1}_pri.bed
		factor=$(wc -l ${1}_pri.bed|awk '{print 1000000/$1}')
		genomeCoverageBed -scale $factor -i ${1}_pri.bed -g $chromsize -bg > ${1}.bdg
		bedSort ${1}.bdg ${1}.bdg1
		bedGraphToBigWig ${1}.bdg1 $chromsize ${1}.bw
		mv ${1}.bdg1 ${1}.bdg
		len=$(gunzip -c $3 |head -n2|tail -n1|awk '{print length($0)}')
		awk -v l=$len -v OFS="\t" '{print $1,$2,$2+l,$4,$5,$6"\n"$1,$3-l,$3,$5,$6}'  ${1}_pri.bed > ${1}_se.bed
	else
		# sam2bam+sort
		samtools view -b -@ $threads -o ${1}.bam ${1}.sam 
		samtools sort -@ $threads -o ${1}_srt.bam ${1}.bam
		# download picard.jar for PE duplicates removal
		wget $picard_url
		# mark duplicates
		java -jar picard.jar MarkDuplicates -I ${1}_srt.bam -O ${1}_mkdup.bam -M ./logs/${1}_dup.log --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT
		echo 'flagstat after mkdup:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_mkdup.bam >> ./logs/${1}_align.log
		# filter our unmapped/failedQC/unpaired/duplicates/secondary alignments
		samtools view -@ $threads -f 2 -F 1804 -b -o ${1}_filtered.bam ${1}_mkdup.bam
		echo >> ./logs/${1}_align.log
		echo 'flagstat after filter:' >> ./logs/${1}_align.log
		samtools flagstat -@ $threads ${1}_filtered.bam >> ./logs/${1}_align.log
		# sort bam by query name for bedpe 
		samtools sort -n -@ $threads -o ${1}.bam2 ${1}_filtered.bam
		bamCoverage --binSize 10 -p $threads --normalizeUsing CPM -b ${1}_filtered.bam -o ${1}.bw
		# bam2bedpe
		bamToBed -bedpe -i ${1}.bam2 > ${1}.bedpe
		bamToBed -i ${1}_filtered.bam > ${1}_se.bed
		# bedpe to standard PE bed for macs2 peak calling (-f BEDPE)
		cut -f1,2,6 ${1}.bedpe > ${1}_pe.bed
		# clean
		rm ${1}_srt.bam ${1}.bam2 ${1}.bedpe picard.jar 
		rm ${1}.bam ${1}.sam
	fi
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "g:x:s:r:i:t:cp:nh" arg
do
	case $arg in
		g) if [ $OPTARG = "hg19" ]; then
			bw2index=$bw2index_hg19
			chromapindex=$chromapindex_hg19
			chromapref=$chromapref_hg19
			chromsize=$chromszize_hg19
		elif [ $OPTARG = "hg38" ]; then
			bw2index=$bw2index_hg38
			chromapindex=$chromapindex_hg38
			chromapref=$chromapref_hg38
			chromsize=$chromszize_hg38
		elif [ $OPTARG = "mm10" ]; then
			bw2index=$bw2index_mm10
			chromapindex=$chromapindex_mm10
			chromapref=$chromapref_mm10
			chromsize=$chromszize_mm10
		else
			echo "Only support hg38, hg19 or mm10, or pass your custom genome build"
			exit 1
		fi;;
		# Bowtie2 index PATH
		x) bw2index=$OPTARG;;
		z) chromsize=$OPTARG;;
		r) chromapref=$OPTARG;;
		i) chromapindex=$OPTARG;;
		t) threads=$OPTARG;;
		c) aln='chromap';;
		p) prefix=$OPTARG;;
		# Nextera adapters
		n) aA='CTGTCTCTTATACACATCT'
		   gG='AGATGTGTATAAGAGACAG';;
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
	prefix=${1%_R1*}
fi

# main
main(){
	if [ ! -d logs ];then 
		mkdir logs
	fi

	if [ ! -d fastqc ];then 
		mkdir fastqc
	fi 
	
	QC_mapping $aln $prefix $1 $2

	sam_bam_bed $prefix $aln $1

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