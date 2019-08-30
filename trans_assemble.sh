#!/bin/bash
#####################################
# Usage:	$READS1=reads1/file/to/path	#
#			$READS2=reads2/file/to/path	#
#			${NAME}=output_file_prefix	#
#			$4=fr|rf|un(null)		#
#####################################
# check programs
requires=("cutadapt" "python3" "hisat2" "fastqc" "stringtie" "STAR")
for i in ${requires[@]};do
	cmd="which "$i" &>/dev/null || { echo \"$i not found\"; exit 1; }"
	eval $cmd 
done

# path to reads
READS1=$1
READS2=$2
NAME=$3
gtf=/genome/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
#gtf=/home/quanyi/app/PLAR/hg19_ensembl.gtf
index=/genome/hg19/HISAT2index/genome_tran
#STAR_idx=/home/quanyi/genome/hg19/STARindex
threads=16

# deternmine the strand direction
if [ $4 = 'fr' ];then
	strand1='--rna-strandness FR'
	strand2='--fr'
	strand3=''
elif [ $4 = 'rf' ];then
	strand1='--rna-strandness RF'
	strand2='--rf'
	strand3=''
else
	strand1='--dta-cufflinks'
	strand2=''
	strand3='--outSAMstrandField intronMotif'
fi

if [ ! -d logs ];then 
mkdir logs
fi

# fastqc control
fastqc -f fastq -t $threads -o logs $READS1 $READS2

# cutadapt--trim adaptors Trueseq index
# python3 version required for -j
cutadapt -m 30 -j $threads -a AGATCGGAAGAGC -A AGATCGGAAGAGC -g GCTCTTCCGATCT -G GCTCTTCCGATCT -o ${NAME}_R1_trimmed.gz -p ${NAME}_R2_trimmed.gz $READS1 $READS2 > ./logs/${NAME}_cutadapt.log

# HISAT2--mapping
hisat2 $strand1 -p $threads -x $index -1 ${NAME}_R1_trimmed.gz -2 ${NAME}_R2_trimmed.gz -S ${NAME}.sam 

samtools view -b -@ $threads ${NAME}.sam -o ${NAME}.bam
samtools sort -@ $threads ${NAME}.bam -o ${NAME}_srt.bam

# STAR--mapping
#STAR --genomeDir /home/quanyi/genome/hg19/STARindex --runThreadN $threads $strand3 --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --outFileNamePrefix ./ --readFilesIn ${NAME}_R1_trimmed.gz ${NAME}_R2_trimmed.gz --readFilesCommand gunzip -c

#stringtie for transcripts assembly
stringtie -G $gtf -l ${NAME} -o ${NAME}.gtf $strand2 ${NAME}_srt.bam
#Aligned.sortedByCoord.out.bam #-B

#clean
rm ${NAME}.bam ${NAME}.sam
mv ${NAME}_srt.bam ${NAME}.bam

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################