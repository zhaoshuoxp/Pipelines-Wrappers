#!/bin/bash
#####################################
# Usage: rmdup_rdm.sh <in bam> <out bam>
#####################################
echo "remove duplicates RANDOMLY with picard"
java -jar /home/quanyi/app/picard.jar MarkDuplicates INPUT=$1 OUTPUT=$2 METRICS_FILE=$(basename $2 .bam)_dup.log REMOVE_DUPLICATES=true DUPLICATE_SCORING_STRATEGY=RANDOM

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################