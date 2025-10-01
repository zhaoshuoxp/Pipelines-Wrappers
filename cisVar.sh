#!/bin/bash
#####################################
# Usage: cisVar.sh <read depth> <individual file>
#####################################
#Check
which cisVar &>/dev/null || { echo "cisVar not found! Install: sudo pip3 install cisVar"; exit 1; }
which samtools &>/dev/null || { echo "samtools not found!"; exit 1; }

#Setup
BAM=$1
DEP=$2
IND=$3
NUM=$(wc -l $IND|cut -f 1 -d ' ')
vcf=/home/quanyi/data/genotypes/72vcf_snp_only/chr*_snps_only.vcf.gz
hg19=/home/quanyi/genome/hg19/GRCh37.p13.genome.fa
prefix=$(basename $BAM .bam)

#Main
cisVar prep -F $prefix -r $DEP -i $IND --chrom-format chr $vcf

cisVar mpileup -F $prefix -r $DEP -f $hg19 -B $BAM 

cisVar post -F $prefix -r $DEP

cisVar geno -F $prefix -r $DEP -i $IND

cisVar qtls -F $prefix -r $DEP -n $NUM

cisVar tidy -F $prefix -r $DEP -b ${prefix}.locations.bed.gz

#Plot pre-post QTLs 
#Need gplots package
plot_fun.R ${prefix}.${DEP}.total.txt

#Clean
rm ${prefix}.genotypes.txt.gz.*.gz.parse.log
################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################