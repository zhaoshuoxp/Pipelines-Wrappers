#!/bin/bash
#####################################
# Usage:                            #
# Manual:                           #
#####################################
which bedToBigBed &>/dev/null || { echo "bedToBigBed not found! Download http://hgdownload.soe.ucsc.edu/admin/exe"; exit 1; }
which sortBed &>/dev/null || { echo "bedtools not found!"; exit 1; }

name=$(basename $1 .txt)
# resolution
res=$2/2
# set q value as color depth
score='$7'
# set  value as color depth if $3="p"
if [ $3 ];then
	if [ $3="p" ];then
		score='$6'
	fi
fi

# convert FitHiC output to interact BED5+13 format
awk -v OFS="\t" '{print $1,$2-'$res',$2+'$res',".",$5,-log('$score'),".", "0",$3,$4-'$res',$4+'$res',".",".", $1,$2-'$res',$2+'$res',".","."}' $1 > $name.temp
# get max -log(q_value), in case of q_value=0
#max_q=$(awk 'BEGIN {max = 0} {if ($6+0> max) max=$6} END {print max}' $name.temp)
# replace "inf" to max -log(q_value)
#awk -v OFS="\t" '{if ($6=="inf")$6='$max_q'}1' $name.temp > $name.temp2
awk -v OFS="\t" '{if ($5>1000)$5=1000}1' $name.temp > $name.temp2

# sort
sortBed -i $name.temp2 > $name.temp3
# get interact.as and chrom size files
wget https://genome.ucsc.edu/goldenPath/help/examples/interact/interact.as
wget http://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes
# convert interact to biginteract
bedToBigBed -as=interact.as -type=bed5+13 $name.temp3 hg19.chrom.sizes $name.bb
   
# clean
rm $name.temp $name.temp2 $name.temp3 interact.as hg19.chrom.sizes 

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################