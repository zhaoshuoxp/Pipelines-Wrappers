#!/bin/bash
#####################################
# Usage: fithic.sh Matrix_raw_contacts BED_bins_coord  bias_IC_normal sample_label      #
# Manual:                           #
#####################################
matrix=$1
abs_bed=$2
bias=$3
sample=$4

# Step 1: covnert HiC-Pro matrix to bed matrix
# hg19 10kb run all chromosomes together 
# pre-process for maximum bin size of 2mb - this will make making bed file and fit-hic pre-processing much easier 
# at 10kb resolution, 2mb = 200 bins, filter for bin2-bin1 < 201
awk '($2-$1)<201' $matrix > $sample.max2m.matrix

# Biases file by hicpro2fithic.py
hicpro2fithic.py -i $matrix -b $abs_bed -s $bias 

run_Mat2Bed.py -c all -k $abs_bed -m $sample.max2m.matrix

# Input 1: Find sum of each group entry in an array: 
# chr1	fragmentMid1	chr2	fragmentMid2	contactCount
awk 'BEGIN{OFS="\t"} {print $1,int(($3+$2)/2),$4,int(($5+$6)/2),$7}' $sample.max2m.matrix.bed > temp_contact.bed

# filter for contact DISTANCE (ie 20kb min, 2mb max)
awk '($4-$2)>20000 && ($4-$2)<2000000' temp_contact.bed > $sample.contact_counts.bed

# Input 2: Filter for contact distances, output fields: 
# chr	extraField	fragmentMid	marginalizedContactCount	mappable
cut -f1,2,5 $sample.contact_counts.bed > temp_frag_left.bed
cut -f3,4,5 $sample.contact_counts.bed > temp_frag_right.bed

cat temp_frag_left.bed temp_frag_right.bed | awk -F"\t" '{a[$1"_"$2]+=$3;} END {for(i in a)  {split(i,b,"_"); print b[1] "\t" 0 "\t" b[2] "\t" a[i] "\t" 1;}}' > $sample.fragment_counts.bed

# Sorting the fragments file is CRITICAL for multi-chromosome inputs
sort -k1,1 -k3,3g $sample.fragment_counts.bed > $sample.fragment_sort_counts.bed

# Make the output directory and run FitHiC
mkdir fithic_out

fithic -f $sample.fragment_sort_counts.bed -i $sample.contact_counts.bed -v -o fithic_out -l $sample.max2m -L 20000 -U 2000000 -t fithic.biases.gz

# q value 0.05 cut-off
zcat *pass2.significances.txt.gz | awk '$7<0.05' > $sample.fithic_q.txt

# Convert to WashU format
awk -v OFS="\t" '{print $1":"$2-5000"-"$2+5000, $3":"$4-5000"-"$4+5000, -log($6)}' $sample.fithic_q.txt > $sample.fithic_q_WashU.txt 

# Clean
rm $sample.max2m.matrix $sample.max2m.matrix.bed temp_contact.bed temp_frag_left.bed temp_frag_right.bed $sample.fragment_counts.bed fithic.fragmentMappability.gz fithic.interactionCounts.gz



################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
#	  Modified from Max Mumbach     #
################ END ################