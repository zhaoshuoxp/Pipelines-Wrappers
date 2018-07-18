#!/bin/bash
#####################################
# Usage:                            #
# Manual:                           #
#####################################
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
# convert FitHiC output to WashU format
awk -v OFS="\t" '{print $1":"$2-'$res'"-"$2+'$res', $3":"$4-'$res'"-"$4+'$res', -log('$score')}' $1 > $1.temp
# get max -log(q_value), in case q_value=0
max_q=$(awk 'BEGIN {max = 0} {if ($3+0> max) max=$3} END {print max}' $1.temp)
# replace "inf" to max -log(q_value)
awk -v OFS="\t" '{if ($3=="inf"){print $1,$2,'$max_q'}else{print $0}}' $1.temp > $(basename $1 .txt)_WashU.txt

rm $1.temp

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################