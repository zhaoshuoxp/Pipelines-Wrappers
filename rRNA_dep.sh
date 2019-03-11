#!/bin/bash
#####################################
# Usage: rRNA_dep.sh <reads1> <reads2> <prefix of output>
#####################################
#rRNA genes BWA index generation
hisat2_rRNA_idx(){
	wd=$(pwd)
	cd /home/quanyi/genome/hg38/HISAT2index/
	grep rRNA $gtf > rRNA.gtf
	gtf_to_fasta rRNA.gtf $genome rRNA.fa
	hisat2-build -p 24 rRNA.fa rRNA
	cd $wd
}

bwa_rRNA_idx(){
	wd=$(pwd)
	cd /home/quanyi/genome/hg38/BWAindex/
	grep rRNA $gtf > rRNA.gtf
	gtf_to_fasta rRNA.gtf $genome rRNA.fa
	bwa index -p rRNA rRNA.fa
	cd $wd
}

# set up
threads=24
genome=/home/quanyi/genome/hg19/GRCh37.p13.genome.fa
gtf=/home/quanyi/genome/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf

if [ ! -f /home/quanyi/genome/hg19/BWAindex/rRNA.fa ]
then
	bwa_rRNA_idx
fi
rRNA_idx=/home/quanyi/genome/hg19/BWAindex/rRNA

# mapping to rRNA index
bwa mem -t 24 $rRNA_idx $1 $2 |samtools view -b -@ 24 -o $3_rRNA.bam

samtools sort -@ 24 -o $3_rRNA_srt.bam $3_rRNA.bam

# count contamination
echo "$3 rRNA contamination(mapped):" > $3_rRNA.log
samtools flagstat -@ 24 $3_rRNA_srt.bam >> $3_rRNA.log

# get unmapped
samtools view -@ 24 -b -f 4 -o $3_umapped.bam $3_rRNA_srt.bam

# get unmapped fastq
samtools fastq -c 6 -@ 24 -1 $3_dep_R1_fq.gz -2 $3_dep_R2_fq.gz $3_umapped.bam

# clean
rm $3_rRNA.bam $3_rRNA_srt.bam $3_umapped.bam

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################