#!/bin/bash
#####################################
# Usage:	RNAseq.sh #
#####################################
# check programs
which cutadapt &>/dev/null || { echo "cutadapt not found!"; exit 1; }
which STAR &>/dev/null || { echo "STAR not found!"; exit 1; }
which fastqc &>/dev/null || { echo "fastqc not found!"; exit 1; }
which featureCounts &>/dev/null || { echo "fastqc not found!"; exit 1; }

gtf=/home/quanyi/genome/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
STAR_idx=/home/quanyi/genome/hg19/STARindex
threads=16

if [ ! -d logs ];then 
mkdir logs
fi

STAR_map(){
	# fastqc control
	fastqc -f fastq -t $threads -o logs $1 
	fastqc -f fastq -t $threads -o logs $2

	# cutadapt--trim adaptors Trueseq index
	# python3 version required for -j
	cutadapt -m 30 -j $threads -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -o $3_R1_trimmed.gz -p $3_R2_trimmed.gz $1 $2 > ./logs/$3_cutadapt.log

	# STAR--mapping
	STAR --genomeDir $STAR_idx --runThreadN $threads --outSAMstrandField intronMotif --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outFileNamePrefix $3 --readFilesIn $1 $2 --readFilesCommand gunzip -c

}


files=(*fasq.gz)
for (( i=0; i<${#files[@]} ; i+=2 ))
do
	prefix=${files[i]%_R1*}
	STAR_map ${files[i]} ${files[i+1]} $prefix
	mv ${prefix}Aligned.sortedByCoord.out.bam ${prefix}.bam
done 

bams=$(ls *.bam)

featureCounts -a $gtf -g gene_name -T $threads -p -t exon -g gene_name -o count.txt $bams

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################