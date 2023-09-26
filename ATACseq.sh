#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
requires=("cutadapt" "python3" "bowtie2" "fastqc" "samtools" "macs2" "bedtools" "bamCoverage" "chromap")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

#### DEFAULT CONFIGURATION ###
# default Paired-end mod
mod='pe'
aln='bowtie2'
# default Nextera adapters
aA='CTGTCTCTTATACACATCT'
gG='AGATGTGTATAAGAGACAG'
# default 1 core to run
threads=1
# genome build url
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
bklt_url_hg19='https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz'
bklt_url_hg38='https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz'
bklt_url_mm10='https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/mm10-blacklist.v2.bed.gz'
# Picard
picard_url='https://github.com/broadinstitute/picard/releases/download/3.1.0/picard.jar'

# help message
help(){
	cat <<-EOF
  Usage: ATAC.sh <options> <reads1>|<reads2> 

  ### INPUT: Paired-end fastq files ###
  This script will QC fastq files and align reads to reference genome build with Bowtie2 or chromap, depending on the species passed by -g or the index and other required files passed by -i, -b and -c, 
  convert alignments to filtered BAM/BED and bigwig,
  then call peaks with MACS2 in BEDPE mode after Tn5 shifting.
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/bowtie2/samtools/bedtools/deeptools/macs2>=2.1.1 required ###

  Options:
    -g [str] Genome build selection <hg38|hg19|mm10>
    -x [str] Custom bowtie2 index PATH
    -b [str] Custom blacklist PATH
    -m [str] Genome size abbr supported by MACS2
    -p [str] Prefix of output
    -t [int] Threads (1 default)
    -s Single-end mod (DO NOT recommend, Paired-end default)
    -c Using chromap to process FASTQ instead of canonical bowtie2
    -i [str] Custom chromap genome index (only valid with -c option)
    -r [str] Custom chromap genome reference (only valid with -c option)
    -z [str] Custom chromosome size table
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
		bowtie2 -X 2000 --local --mm -p $threads -x $bw2index -U ${2}_trimmed.fastq.gz -S ${2}.sam
		echo 'Bowtie2 mapping summary:' > ./logs/${2}_align.log
		tail -n 15 nohup.out >> ./logs/${2}_align.log
	else
		# paired-end CMD
		# FastQC
		fastqc -f fastq -t $threads -o fastqc $3 $4
		# TruSeq adapter trimming
		cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${2}_trimmed_R1.fastq.gz -p ${2}_trimmed_R2.fastq.gz $3 $4 > ./logs/${2}_cutadapt.log
		# Bowtie2 align
		bowtie2 -X 2000 --local --mm -p $threads -x $bw2index -1 ${2}_trimmed_R1.fastq.gz -2 ${2}_trimmed_R2.fastq.gz -S ${2}.sam
		echo 'Bowtie2 mapping summary:' > ./logs/${2}_align.log
		tail -n 15 nohup.out >> ./logs/${2}_align.log
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
		wget $picard_url
		# mark duplicates
		java -jar picard.jar MarkDuplicates -I ${1}_srt.bam -O ${1}_mkdup.bam -M ./logs/${1}_dup.log --REMOVE_DUPLICATES false --VALIDATION_STRINGENCY SILENT
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
		bamToBed -i ${1}_filtered.bam > ${1}.bed
		# clean
		rm ${1}_srt.bam ${1}.bam2 picard.jar 
	fi
	rm ${1}.bam ${1}.sam ${1}_chrM.bam 
	# convert filtered BAM to CPM normalized bigWig with deeptools
	bamCoverage --binSize 10 -p $threads --normalizeUsing CPM -b ${1}_filtered.bam -o ${1}.bw
}

# Peak calling with MACS2 >v2.1.1
peak_calling(){
	if [ $1 = 'se' ];then
		# Tn5 shift in SE mode
		awk -F $'\t' 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4}else if($6=="-"){$3=$3-5} print $0}' ${2}_se.bed > ${2}_shift_se.bed
		mv ${2}_shift_se.bed  ${2}_se.bed
		# broad peak calling
		cd macs2
		macs2 callpeak -t ../${2}_se.bed -g $sp -n ${2} -f BED --keep-dup all --broad --nomodel --shift -37 --extsize 73
		# Blacklist filter 
		intersectBed -v -a ${2}_peaks.broadPeak -b $blkt_file > ${2}_broad_filtered.bed
		cd ..
	else
		# Tn5 shift in PE mode
		awk -v OFS="\t" '{if($9=="+"){print $1,$2+4,$6+4}else if($9=="-"){if($2>=5){print $1,$2-5,$6-5}else if($6>5){print $1,0,$6-5}}}' ${2}.bedpe > ${2}_pe.bed
		awk -F $'\t' 'BEGIN{OFS=FS}{if($6=="+"){$2=$2+4}else if($6=="-"){$3=$3-5} print $0}' ${2}.bed > ${2}_se.bed
		# broad peak calling
		cd macs2
		echo "MACS2 version >= 2.1.1 required!"
		macs2 callpeak -t ../${2}_pe.bed -g $sp -n ${2} -f BEDPE --keep-dup all --broad 
		# Blacklist filter 
		intersectBed -v -a ${2}_peaks.broadPeak -b $blkt_file > ${2}_broad_filtered.bed
		cd ..
	fi
	rm ${2}.bedpe ${2}.bed 
}

chromap_total(){
	if [ $2 = 'se' ];then
		chromap --preset atac -r $chromapref -x $chromapindex -t $threads -1 $3 -o ${1}_se.bed
		echo 'chromap mapping summary:' > ./logs/${1}_align.log
		tail -n 14 nohup.out >> ./logs/${1}_align.log
		awk 'substr($1,1,3)=="chr"' ${1}_se.bed > ${1}_pri.bed
	else
		chromap --preset atac -r $chromapref -x $chromapindex -t $threads -1 $3 -2 $4 -o ${1}_pe.bed
		echo 'chromap mapping summary:' > ./logs/${1}_align.log
		tail -n 14 nohup.out >> ./logs/${1}_align.log
		awk 'substr($1,1,3)=="chr"' ${1}_pe.bed > ${1}_pri.bed
		len=$(gunzip -c $3 |head -n2|tail -n1|awk '{print length($0)}')
		awk -v l=$len -v OFS="\t" '{print $1,$2,$2+l,$4,$5,$6"\n"$1,$3-l,$3,$4,$5,$6}'  ${1}_pri.bed > ${1}_se.bed
	fi
	factor=$(wc -l ${1}_pri.bed|awk '{print 1000000/$1}')
	genomeCoverageBed -scale $factor -i ${1}_pri.bed -g $chromsize -bg > ${1}.bdg
	bedSort ${1}.bdg ${1}.bdg1
	bedGraphToBigWig ${1}.bdg1 $chromsize ${1}.bw
	mv ${1}.bdg1 ${1}.bdg

	cd macs2
	echo "MACS2 version >= 2.1.1 required!"
	if [ $2 = 'se' ];then
		mv ../${1}_pri.bed ../${1}_se.bed
		macs2 callpeak -t ../${1}_se.bed -g $sp -n ${1} -f BED --keep-dup all --broad
	else
		mv ../${1}_pri.bed ../${1}_pe.bed 
		macs2 callpeak -t ../${1}_pe.bed -g $sp -n ${1} -f BEDPE --keep-dup all --broad 
	fi
	# Blacklist filter 
	intersectBed -v -a ${1}_peaks.broadPeak -b $blkt_file > ${1}_broad_filtered.bed
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "g:x:b:m:t:sp:z:i:r:ch" arg
do
	case $arg in
		g) if [ $OPTARG = "hg19" ]; then
			bw2index=$bw2index_hg19
			chromapindex=$chromapindex_hg19
			chromapref=$chromapref_hg19
			chromsize=$chromszize_hg19
			curl -s $bklt_url_hg19 | gunzip -c > bklt
			sp='hs'
		   elif [ $OPTARG = "hg38" ]; then
			bw2index=$bw2index_hg38
			chromapindex=$chromapindex_hg38
			chromapref=$chromapref_hg38
			chromsize=$chromszize_hg38
			curl -s $bklt_url_hg38 | gunzip -c > bklt
			sp='hs'
		   elif [ $OPTARG = "mm10" ]; then
			bw2index=$bw2index_mm10
			chromapindex=$chromapindex_mm10
			chromapref=$chromapref_mm10
			chromsize=$chromszize_mm10
			curl -s $bklt_url_mm10 | gunzip -c > bklt
			sp='mm'
		   else
			echo "Only support hg38, hg19 or mm10, or pass your custom genome build"
			exit 1
		   fi
		   blkt_file=$(readlink -f bklt);;
		# Bowtie2 index PATH
		x) bw2index=$OPTARG;;
		b) blkt_file=$OPTARG;;
		m) sp=$OPTARG;;
		t) threads=$OPTARG;;
		# single-end mod
		s) mod='se';;
		p) prefix=$OPTARG;;
		z) chromsize=$OPTARG;;
		r) chromapref=$OPTARG;;
		i) chromapindex=$OPTARG;;
		c) aln='chromap';;
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
	
	if [ ! -d macs2 ];then 
		mkdir macs2
	fi 
	
	if [ $aln = 'chromap' ];then
		chromap_total $prefix $mod $1 $2
	else
		QC_mapping $mod $prefix $1 $2

		sam_bam_bed $prefix $mod

		peak_calling $mod $prefix
	fi

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