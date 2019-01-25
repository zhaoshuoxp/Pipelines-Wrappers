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
# convert FitHiC output to WashU longrange format
awk -v OFS="," '{print $1"\t"$2-'$res'"\t"$2+'$res'"\t"$3":"$4-'$res'"-"$4+'$res',-log('$score')"\n"$3"\t"$4-'$res'"\t"$4+'$res'"\t"$1":"$2-'$res'"-"$2+'$res',-log('$score')}' $1 > $1.temp
# get max -log(q_value), in case of q_value=0
max_q=$(awk -F"," 'BEGIN {max = 0} {if ($2+0> max) max=$2} END {print max}' $1.temp)
# replace "inf" to max -log(q_value)
awk -v OFS="\t" -F"," '{if ($2=="inf"){print $1,'$max_q'}else{print $0}}' $1.temp > $(basename $1 .txt).longrange
# sort file
sortBed -i $(basename $1 .txt).longrange > $(basename $1 .txt).lr
# bgzip compress and index
bgzip $(basename $1 .txt).lr
tabix -p bed $(basename $1 .txt).lr.gz
# clean
rm $1.temp $(basename $1 .txt).longrange

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################