#!/bin/bash
#####################################
# Usage:                            #
# Manual:                           #
#####################################
#rRNA genes BWA index generation
bwa_rRNA_idx(){
	genome=/home/quanyi/genome/hg19/GRCh37.p13.genome.fa
	gtf=/home/quanyi/genome/hg19/gencode.v19.chr_patch_hapl_scaff.annotation.gtf
	grep rRNA $gtf > rRNA.gtf
	gtf_to_fasta rRNA.gtf $genome rRNA.fa
	bwa index rRNA.fa
}

# set up
thread=32
rRNA_idx=/home/quanyi/genome/hg19/BWArRNA/rRNA.fa

# mapping to rRNA index
bwa mem -t $thread $rRNA_idx $1 $2 |samtools view -b -@ $thread -o $3_rRNA.bam

# sort BAM
samtools sort -@ $thread -o $3_rRNA_srt.bam $3_rRNA.bam

# mapped rate = rRNA contamination
echo "$3 rRNA contamination:" > $3_rRNA.log
samtools flagstat -@ $thread $3_rRNA_srt.bam >> $3_rRNA.log



################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################