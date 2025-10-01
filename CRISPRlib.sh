#!/bin/bash
#####################################
requires=("cutadapt" "python3" "bowtie" "samtools")
for i in ${requires[@]};do
    which $i &>/dev/null || { echo $i not found; exit 1; }
done


help(){
    cat <<-EOF
    Usage: CRISPRlib.sh <options> <reads_clean.fq.gz> 

    ### INPUT: fastq files ###
    This script will trim the input fastq to 20nt after the given sequence with cutadapt, and align the trimmed reads to the reference library build with Bowtie1, depending on the library selection passed by -l or the index and adapter sequence passed by -i and -a,
    then statisticize each sequence's frequency, and all results will be store in current (./) directory.
    ### python3/cutadapt/bowtie1/samtools required ###

    Options:
    -l [str] library selection <110066|160129|162256>
    -i [str] Custom bowtie index PATH
    -a [str] Custom adapter sequence
    -p [str] Prefix of output
    -n [int] Threads (1 default)
    -h Print this help message
EOF
    exit 0
}

main(){
    cutadapt -g $adpt -j $threads -l 20 -m 19 -o ${pre}tr.fq.gz $1 > ${pre}log
    
    bowtie -x $idx -n 0 -p $threads --no-unal -l 20 ${pre}tr.fq.gz -S ${pre}sam 2>&1|tee -a ${pre}log
    
    samtools view -@ $threads -o ${pre}bam ${pre}sam
    samtools sort -@ $threads -o ${pre}srt.bam ${pre}bam
    samtools index -@ $threads ${pre}srt.bam
    
    samtools idxstats  ${pre}srt.bam |awk '$3+0>0' |awk '{print $1"\t"$3}' > ${pre}counts.tsv
    
    awk 'substr($1,1,1)!="@"' ${pre}sam|awk '{print $3"\t"$10}'|sort |uniq -c |awk '{print $2"\t"$3"\t"$1}' > ${pre}table.tsv
}

# no ARGs error
if [ $# -lt 1 ];then
    help
    exit 1
fi


while getopts "l:i:a:n:p:h" arg
do
    case $arg in
        l) if [ $OPTARG = "110066" ]; then
            idx='/mnt/date3/Project/zhaoqy/genome/CRISPR/hs_metabolism_110066/metabolism'
            adpt='TTTCTAGCTCTAAAAC'
        elif [ $OPTARG = "160129" ]; then
            idx='/mnt/date3/Project/zhaoqy/genome/CRISPR/mm_metabolism_160129/sgRNA_lib'
            adpt='TGTTTCCAGCATAGCTCTTAAAC'
        elif [ $OPTARG = "162256" ]; then
	    idx='/mnt/date3/Project/zhaoqy/genome/CRISPR/hs_epigentics_162256/lib'
	    adpt='ATAGCTCTTAAAC'
	else
            echo "Only support 110066, 160129, 162256 libraries, or pass your own bowtie1 index and adapter sequences"
            exit 1
        fi;;
        i) idx=$OPTARG;;
        a) adpt=$OPTARG;;
        n) threads=$OPTARG;;
        p) pre=$OPTARG;;
        h) help ;;
        ?) help
            exit 1;;
    esac
done

# shift ARGs to reads
shift $(($OPTIND - 1))
# get prefix of output
if [ -z $pre ];then
    echo "No -p <prefix> given, use file name as prefix"
    pre=${1/clean.fq.gz/}
fi

if [ -z $threads ];then
    echo "Using 12 threads as default"
    threads=12
fi

main $1

# check running status
if [ $? -ne 0 ]; then
    help
    exit 1
else
    echo "Run succeed"
fi

