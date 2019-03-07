#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
which cutadapt &>/dev/null || { echo "cutadapt not found, install by pip3!"; exit 1; }
which python3 &>/dev/null || { echo "python3 not found, install python3!"; exit 1; }

# help message
help(){
	cat <<-EOF
  adapt_trim.sh <options> <reads1>|..<reads2> 
  Trim adapter sequences from fastq files with cutadapt@python3.
  Trimmed and gzipped fastq will be store in current (./) directory.
  Options:
    -p Prefix of output
    -t Threads (1 default)
    -s Single-end mod (Paired-end default)
    -n Nextera adapters (Truseq default)
    -h Print this help message
EOF
	exit 0
}

# cutadapt CMD
adapt_trim(){
	if [ $1 = 'se' ];then
		# single-end CMD
		cutadapt -m 30 -j $2 -a $3 -g $4 -o ${5}_trimmed.fastq.gz $6 > ${5}_cutadapt.log
	else
		# paired-end CMD
		cutadapt -m 30 -j $2 -a $3 -A $3 -g $4 -G $4 -o ${5}_trimmed_R1.fastq.gz -p ${5}_trimmed_R2.fastq.gz $6 $7 > ${5}_cutadapt.log
	fi
}

# no ARGs error
if [ $# -lt 1 ]
then
	help
	exit 1
fi

# default Paired-end mod
mod='pe'
# default TrueSeq adapters
aA='AGATCGGAAGAGC'
gG='GCTCTTCCGATCT'
# default 1 core to run
threads=1

while getopts "snt:hp:" arg
do
	case $arg in
		t) threads=$OPTARG;;
		# single-end mod
		s) mod='se';;
		# Nextera adapters
		n) aA='CTGTCTCTTATACACATCT'
		   gG='AGATGTGTATAAGAGACAG';;
		p) prefix=$OPTARG;;
		h) help ;;
		?) help
			exit 1;;
	esac
done

# shift ARGs to reads
shift $(($OPTIND - 1))
# get prefix of output
if [ -z $prefix ];then
	echo "No -p <prefix> given, use file name as prefix"
	if [ $mod = 'se' ];then
		prefix=${1%.*}
	else
		prefix=${1%_R1*}
	fi
fi

# main
adapt_trim $mod $threads $aA $gG $prefix $1 $2

# check running status
if [ $? -ne 0 ]; then
	help
	exit 1
else
	echo "Run succeed"
fi

################ END ################
#          Created by Aone          #
#     quanyi.zhao@stanford.edu      #
################ END ################