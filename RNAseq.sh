#!/bin/bash

# check programs
# multi-core support requires cutadapt installed and run by python3
requires=("cutadapt" "python3" "STAR" "fastqc" "featureCounts")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

#### DEFAULT CONFIGURATION ####
# default TruSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1
# STAR index
STAR_idx_hg='/nfs/baldar/quanyiz/genome/hg38/STARindex'
STAR_idx_mm='/nfs/baldar/quanyiz/genome/mm10/STARindex'
# GTF ref
gtf_hg='/nfs/baldar/quanyiz/genome/hg38/gencode.v44.chr_patch_hapl_scaff.annotation.gtf'
gtf_mm='/nfs/baldar/quanyiz/genome/mm10/gencode.vM25.chr_patch_hapl_scaff.annotation.gtf'

# help message
help(){
	cat <<-EOF
  Usage: RNAseq.sh <options> -m meta.txt </PATH/to/fastq/>

  ### INPUT: Paired-end fastq files with _R1/2.fastq.gz suffix, and a text file with meta information  ###
  This script will QC fastq files and align reads to the reference genome and transcriptome with STAR, depending on the species  passed by -s or the index and GTF passed by -i and -g,
  featureCounts and DESeq2 will be used for reads counting and differentially expressed genes test,
  All results will be store in current (./) directory.
  ### indice and GTF have to be the same assembly version ###
  ### python3/cutadapt/fastqc/STAR>=2.7/R/featureCounts/DEseq2 required ###

  Options:
    -m [str] /PATH/to/meta.txt
    -s [str] species <hg|mm> hg=hg38, mm=mm10
    -i [str] Custom STAR index PATH
    -g [str] Custom Reference GTF transcripts PATH
    -t [int] Threads (1 default)
    -p prepare meta.txt
    -n Nextera adapters (Truseq default)
    -h Print this help message

  NOTE:
    1) ### Put fastq files and give the their parent directory (NOT fastq files themselves, script will scan the directory) to the script ###
    2) Sample names in meta.txt must be shown without _R1/2.fastq.gz suffix, The order of samples has to the same as in command: ls -1 for the script to work,
    You may use this script to prepare the meta.txt:
      RNAseq.sh -p </PATH/to/fastq/directory/>
    Then edit meta.txt generated in current directory to add group/condition in the 2nd column,
    Provide this text to the script by <-m meta.txt>.
EOF
	exit 0
}

STAR_map(){
	# fastqc control
	fastqc -f fastq -t $threads -o logs $1 $2
	# cutadapt--trim adaptors Truseq index
	# python3 version required for -j
	cutadapt -m 30 -j $threads -a $aA -A $aA -o ${3}_R1_trimmed.gz -p ${3}_R2_trimmed.gz $1 $2 > ./logs/${3}_cutadapt.log
	#cutadapt -m 30 -j $threads -a $aA -A $aA -g $gG -G $gG -o ${3}_R1_trimmed.gz -p ${3}_R2_trimmed.gz $1 $2 > ./logs/${3}_cutadapt.log
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

while getopts "m:t:s:g:nhi:p" arg
do
	case $arg in
		m) meta=$OPTARG;;
		t) threads=$OPTARG;;
		s) if [ $OPTARG = "hg" ]; then
			gtf=$gtf_hg
			STAR_idx=$STAR_idx_hg
		   elif [ $OPTARG = "mm" ]; then
			gtf=$gtf_mm
			STAR_idx=$STAR_idx_mm
		   else
			echo "Only support hg or mm, or pass your custom STAR index and GTF by -i and -g"
			exit 1
		   fi;;
		g) gtf=$OPTARG;;
		# Nextera adapters
		n) aA='CTGTCTCTTATACACATCT'
		   gG='AGATGTGTATAAGAGACAG';;
		i) STAR_idx=$OPTARG;;
		p) shift $(($OPTIND - 1))
		   echo -e "Sample\tGroup" > meta.txt
		   files=($1/*.f*q.gz)
		   for (( i=0; i<${#files[@]} ; i+=2 ))
				do
					filename=${files[i]##*/}
					prefix=${filename%_*1.f*q.gz}
					echo -e "$prefix\t" >> meta.txt
				done
		   vim meta.txt
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

	files=($1/*.f*q.gz)
	for (( i=0; i<${#files[@]} ; i+=2 ))
	do
		filename=${files[i]##*/}
		prefix=${filename%_*1.f*q.gz}
		STAR_map ${files[i]} ${files[i+1]} $prefix
		rm -r ${prefix}_STARtmp
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
	mv *SJ.out.tab BAM/

	cat >deseq.r<<-EOF
	#!/usr/bin/env Rscript
	library("DESeq2")
	options<-commandArgs(trailingOnly = T)

	data <- read.table(options[1], header=T, quote="\t", check.names=F, skip=1)
	meta <- read.table(options[2], header=T, quote="\t")
	
	counTpm <- function(counts, lengths) {
		rate <- counts / lengths
		rate / sum(rate) * 1e6 }
	
	names(data)[7:ncol(data)] ->sample
	sub(pattern=".bam", replacement="", sample) ->sampleNames
	if (all(meta\$Sample==sampleNames)){

		names(data)[7:ncol(data)] <- sampleNames
		countData <- as.matrix(data[7:ncol(data)])
		rownames(countData) <- data\$Geneid

		# export TPM values for all genes
		counTpm(countData,data\$Length)->tpm
		write.table(tpm, "allgenes_TPM.txt", row.names=T, sep="\t", quote=F)

		database <- data.frame(name=sampleNames, group=meta\$Group)
		dds <- DESeqDataSetFromMatrix(countData, colData=database, design= ~ group)
		dds <- dds[ rowSums(counts(dds)) > 1, ]
		dds <- DESeq(dds)
		res <- results(dds)

		# Only export normalized Exp values
		#res <- results(dds, contrast=c('group',"treat","control"))
		#write.table(res, "DEG.txt")
		resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)
		#resdata <- as.data.frame(counts(dds, normalized=TRUE))
		names(resdata)[names(resdata)=="Row.names"]="Genes"
		write.table(resdata, "allgenes_DESeq2.txt", row.names=F, sep="\t", quote=F)

	}else{
		stop("Sample names in meta.txt don't match featureCount output!")
	}
EOF
	chmod 755 deseq.r
	./deseq.r featureCounts.txt $2
	#rm deseq.r
}

main $1 $meta

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