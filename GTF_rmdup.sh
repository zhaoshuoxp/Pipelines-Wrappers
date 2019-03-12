#!/bin/bash
#####################################

gtf=$1
pre=$(basename $gtf .gtf)

gtfToGenePred $gtf ${pre}.pd
genePredToBed ${pre}.pd ${pre}.bed

num_before=$(wc -l ${pre}.bed|awk '{print $1}')
echo "Total $num_before transcripts before deduplication"

echo "Remove duplicates by random"
sort -k1,1 -k2,3n -k11,12 -u ${pre}.bed > ${pre}_uniq.bed

cut -f4 ${pre}_uniq.bed > uniq.list
num_after=$(wc -l uniq.list|awk '{print $1}')
echo "Total $num_after transcripts after deduplication"

echo "Writing deduplicated GTF..."
awk -F'"' 'NR==FNR{a[$0]}NR>FNR{if ($4 in a) print $0}' uniq.list $gtf > ${pre}_uniq.gtf
echo "Done"

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################