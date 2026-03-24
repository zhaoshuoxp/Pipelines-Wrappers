#!/bin/bash
#####################################

requires=("cutadapt" "python3" "bowtie" "bowtie-build" "samtools")
for i in ${requires[@]};do
    which $i &>/dev/null || { echo $i not found; exit 1; }
done

help(){
    cat <<-EOF
    Usage: CRISPRlib.sh <options> <reads_clean.fq.gz> 

    ### INPUT: fastq files ###
    This script takes an sgRNA list (TSV/CSV), builds a Bowtie1 index,
    trims the input fastq to 20nt after the given adapter sequence with cutadapt, 
    aligns the trimmed reads to the newly built reference library with Bowtie1,
    then statisticizes each sequence's frequency. All results will be stored in current (./) directory.
    ### python3/cutadapt/bowtie1/samtools required ###

    Options:
    -s [file] sgRNA sequence list (TSV or CSV format: ID <tab/comma> sequence) [Required]
    -a [str]  Custom adapter sequence [Required]
    -p [str]  Prefix of output
    -n [int]  Threads (12 default)
    -h        Print this help message
EOF
    exit 0
}


build_index(){
    local sgrna_file=$1
    local prefix=$2
    echo "Converting ${sgrna_file} to FASTA and building Bowtie index..."

    awk -F '[, \t]+' '{gsub(/\r/,"",$2); print ">"$1"\n"$2}' "$sgrna_file" > "${prefix}ref.fa"

    bowtie-build "${prefix}ref.fa" "${prefix}ref_index" --quiet
    echo "Index built successfully."
}

main(){
    local fastq=$1
    
    echo "Running cutadapt..."
    cutadapt -g $adpt -j $threads -l 20 -m 19 -o ${pre}tr.fq.gz "$fastq" > ${pre}log
    
    echo "Running bowtie alignment..."
    bowtie -x "${pre}ref_index" -n 0 -p $threads --no-unal -l 20 ${pre}tr.fq.gz -S ${pre}sam 2>&1 | tee -a ${pre}log
    
    echo "Processing SAM/BAM files..."
    samtools view -@ $threads -b -o ${pre}bam ${pre}sam  # 添加了 -b 确保输出 BAM
    samtools sort -@ $threads -o ${pre}srt.bam ${pre}bam
    samtools index -@ $threads ${pre}srt.bam
    
    echo "Calculating statistics..."
    samtools idxstats ${pre}srt.bam | awk '$3+0>0' | awk '{print $1"\t"$3}' > ${pre}counts.tsv
    
    awk 'substr($1,1,1)!="@"' ${pre}sam | awk '{print $3"\t"$10}' | sort | uniq -c | awk '{print $2"\t"$3"\t"$1}' > ${pre}table.tsv
}

if [ $# -lt 1 ];then
    help
fi

while getopts "s:a:n:p:h" arg
do
    case $arg in
        s) sgrna=$OPTARG;;
        a) adpt=$OPTARG;;
        n) threads=$OPTARG;;
        p) pre=$OPTARG;;
        h) help ;;
        ?) help
            exit 1;;
    esac
done

shift $(($OPTIND - 1))

fastq_input=$1
if [ -z "$fastq_input" ]; then
    echo "Error: Missing input FastQ file."
    help
fi

if [ -z "$sgrna" ] || [ -z "$adpt" ]; then
    echo "Error: Both -s <sgRNA file> and -a <adapter> are required."
    exit 1
fi

if [ -z "$pre" ];then
    echo "No -p <prefix> given, use file name as prefix"
    pre=${fastq_input/clean.fq.gz/}
fi

if [ -z "$threads" ];then
    echo "Using 12 threads as default"
    threads=12
fi

build_index "$sgrna" "$pre"

main "$fastq_input"

if [ $? -ne 0 ]; then
    echo "Run failed"
    exit 1
else
    echo "Run succeed"
fi