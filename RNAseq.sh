#!/bin/bash

# check programs
# multi-core support requires cutadapt installed and run by python3
which cutadapt &>/dev/null || { echo "cutadapt not found, install by pip3!"; exit 1; }
which python3 &>/dev/null || { echo "python3 not found, install python3!"; exit 1; }
which STAR &>/dev/null || { echo "STAR not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which featureCounts &>/dev/null || { echo "fastqc not found!"; exit 1; }

#### DEFAULT CONFIGURATION #### 

# default TruSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1
# STAR index
STAR_idx='/home/quanyi/genome/hg19/STARindex'
# GTF ref
gtf='/home/quanyi/genome/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf'

# help message
help(){
	cat <<-EOF
  Usage: RNAseq.sh <options> -c conditions.txt </PATH/to/fastq/> 

  ### INPUT: Paired-end fastq files with _R1/2.fastq.gz extension, and a text file discribing samples per conditon ###
  This script QC fastq files and align reads to hg19/GRCh37(depends on index and GTF provided) using STAR, 
  featureCounts and DESeq2 will be used for reads count and differntial expresss genes discovery,
  All results will be store in current (./) directory.
  ### python3/cutadapt/fastqc/STAR/R/featureCounts/DEseq2 required ###

  Options:
    -c [str] /PATH/to/conditions.txt
    -i [str] STAR index PATH
    -g [str] Reference GTF transcripts PATH
    -t [int] Threads (1 default)
    -p prepare condition.txt
    -n Nextera adapters (Truseq default)
    -h Print this help message

  NOTE:
    1) ### Put fastq files of each condition all together in a directoy and give this PATH (NOT fastq files) to the script ###
    2) Sample names in conditions.txt must be shown without _R1/2.fastq.gz extension,
 	The order of the samples has to the same as in command: ls -1 for the script to work,
    You may use this script to prepare the conditions.txt:
      RNAseq.sh -p </PATH/contains/fastq> 
    Then edit conditions.txt in current directory by adding condition names in 2nd column,
    Provide this text to the script by <-c conditions.txt>.
EOF
	exit 0
}

STAR_map(){
	# fastqc control
	fastqc -f fastq -t $threads -o logs $1 $2
	# cutadapt--trim adaptors Truseq index
	# python3 version required for -j
	cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${3}_R1_trimmed.gz -p ${3}_R2_trimmed.gz $1 $2 > ./logs/${3}_cutadapt.log
	# set open file limit for STAR BAM sorting 
	ulimit -n 10000
	# STAR--mapping
	STAR --genomeDir $STAR_idx --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $3 --readFilesIn ${3}_R1_trimmed.gz ${3}_R2_trimmed.gz --readFilesCommand gunzip -c
	mv ${3}Log.out logs
	mv ${3}Log.progress.out logs
	mv ${3}Log.final.out logs
	mv ${3}Aligned.sortedByCoord.out.bam ${3}.bam
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "c:t:g:nhi:p" arg
do
	case $arg in
		c) cond=$OPTARG;;
		t) threads=$OPTARG;;
		# single-end mod
		g) gtf=$OPTARG;;
		# Nextera adapters
		n) aA='CTGTCTCTTATACACATCT'
		   gG='AGATGTGTATAAGAGACAG';;
		i) STAR_idx=$OPTARG;;
		p) shift $(($OPTIND - 1))
		   echo -e "sample\tcondition" > conditions.txt
		   files=($1/*fastq.gz)
		   for (( i=0; i<${#files[@]} ; i+=2 ))
				do
					filename=${files[i]##*/}
					prefix=${filename%_R1*}
					echo -e "$prefix\t" >> conditions.txt
				done
		   vim conditions.txt 
		   exit 0;;
		h) help ;;
		?) help
			exit 1;;
	esac
done

# shift ARGs to reads
shift $(($OPTIND - 1))

main(){
	if [ ! -d logs ];then 
		mkdir logs
	fi

	files=($1/*fastq.gz)
	for (( i=0; i<${#files[@]} ; i+=2 ))
	do
		filename=${files[i]##*/}
		prefix=${filename%_R1*}
		STAR_map ${files[i]} ${files[i+1]} $prefix
		bam=${bam}" "${prefix}.bam
	done 
	
	if [ ! -d TRIMMED ];then 
		mkdir TRIMMED
	fi
	mv *_trimmed.gz TRIMMED/
	
	featureCounts -a $gtf -g gene_name -T $threads -p -t exon -o featureCounts.txt $bam
	
	if [ ! -d BAM ];then 
		mkdir BAM
	fi
	mv *.bam BAM/
	mv *.SJ.out.tab BAM/
	
	cat >deseq.r<<-EOF
	#!/usr/bin/env Rscript
	library("DESeq2")
	options<-commandArgs(trailingOnly = T)
	
	data <- read.table(options[1], header=T, quote="\t", check.names=F, skip=1)
	SampCond <- read.table(options[2], header=T, quote="\t")

	names(data)[7:ncol(data)] ->sample
	sub(pattern=".bam", replacement="", sample) ->sampleNames
	if (all(SampCond\$sample==sampleNames)){
		
		names(data)[7:ncol(data)] <- sampleNames
		countData <- as.matrix(data[7:ncol(data)])
		rownames(countData) <- data\$Geneid
		database <- data.frame(name=sampleNames, condition=SampCond\$condition)
		rownames(database) <- sampleNames
		
		dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ condition)
		dds <- dds[ rowSums(counts(dds)) > 1, ]
		dds <- DESeq(dds)
		res <- results(dds)
		#write.table(res, "DEG.txt")
		resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
		names(resdata)[names(resdata)=="Row.names"]="Genes"
		write.table(resdata, "all_genes_exp.txt", row.names=F, sep="\t", quote=F)

	}else{
		stop("Sample names in conditions.txt don't match featureCount output!")	
	}
	EOF
	chmod 755 deseq.r 
	./deseq.r featureCounts.txt $2	
	rm deseq.r
}

main $1 $cond

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