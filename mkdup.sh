#!/bin/bash
#####################################
# Usage:                            #
# Manual:                           #
#####################################

java -jar /home/quanyi/app/picard.jar MarkDuplicates INPUT=$1 OUTPUT=$2 METRICS_FILE=$3 REMOVE_DUPLICATES=false

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################