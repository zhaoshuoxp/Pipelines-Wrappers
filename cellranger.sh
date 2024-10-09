#!/bin/bash

# check dependences
# multi-core support requires cutadapt installed and run by python3
requires=("cellranger-atac" "cellranger-arc" "cellranger" "cellranger7.2")
for i in ${requires[@]};do
	which $i &>/dev/null || { echo $i not found; exit 1; }
done

#### DEFAULT CONFIGURATION ###
# default running mode
mem=200
threads=20
norm='none'
aggr='none'
# cellranger reference path
arc_path='/home/quanyiz/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0/'
rna_path='/home/quanyiz/genome/refdata-gex-mm10-2020-A_tdT/'
cellranger_path='cellranger7.2'
ref_type='--transcriptome'
ref_path=$rna_path
# help message
help(){
	cat <<-EOF
  Usage: cellranger.sh <options> <path/to/fastqs|output/>

  ### INPUT: Paired-end FASTQ files following the cellranger demultiplexed naming conventions ###
	This script runs cellranger automatically for all samples in the specified fastq directory.
	Sample names are inferred from the fastq filenames, i.e <prefix>_S1_L001_R1_001.fastq.gz. 
	Outputs are saved in the current directory with filenames based on the inferred prefix. 
	The script supports 10x genomic RNA (default), ATAC, and multiome based on the specified data type (-m).
  Options:
    -m [str] Data type selection <rna|atac|multiome>, required in both counting and aggregation
    -x [str] Custom cellranger reference PATH
    -t [int] Threads (20 by default)
    -r [int] RAM used (200 gig by default)
    -u Use cellranger8.0.1 in RNA counting, only work with -m rna
    -a Run cellranger "aggr" function instead of "count", if no CSV file provided, a CSV file will be automatically generated by scanning all the sample output folders in the provided path (Only work with ATAC and multiome data)
    -c [str] Custom CSV file for cellranger aggragation, only work with -a
    -n Normalize samples by sequencing depth in aggregation function, only work with -a (NOT performed by default)
    -h Print this help message
EOF
	exit 0
}

# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "m:x:t:r:uac:nh" arg
do
	case $arg in
		m) if [ $OPTARG = "rna" ]; then
			cellranger_path='cellranger7.2'
			ref_type='--transcriptome'
			ref_path=$rna_path
			mod='rna'
		   elif [ $OPTARG = "atac" ]; then
			cellranger_path='cellranger-atac'
			ref_type='--reference'
			ref_path=$arc_path
			mod='atac'
		   elif [ $OPTARG = "multiome" ]; then
			cellranger_path='cellranger-arc'
			ref_type='--reference'
			ref_path=$arc_path
			mod='arc'
		   else
			echo "Only support rna, atac or multiome"
			exit 1
		   fi;;
		#custom ref_oath
		x) ref_path=$OPTARG;;
		t) threads=$OPTARG;;
		r) mem=$OPTARG;;
		# use cellranger8.0.1
		u) cellranger_path='cellranger'
		   ref_type='--transcriptome'
		   ref_path=$rna_path
		   added_par='--create-bam true';;
		a) aggr='aggr';;
		c) csv=$OPTARG;;
		n) norm='depth';;
		h) help ;;
		?) help
			exit 1;;
	esac
done

# shift ARGs to reads
shift $(($OPTIND - 1))

# main
main(){
	if [ -z $mod ];then
		echo "Data type requeired by -m"
		exit 1
	fi
	if [ $aggr = "aggr" ]; then
		if [ -z $csv ];then
			if [ mod = "atac" ];then
				echo "library_id,fragments,cells" > aggr.csv
				ls -F | grep "/$"|cut -f1 -d'/' > id
				ls -1 */outs/fragments.tsv.gz > frag
				ls -1 */outs/singlecell.csv > cell
				paste id frag cell -d',' >>aggr.csv
				csv='aggr.csv'
				rm id cell frag
			elif [ mod = "arc" ]; then
				echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > aggr.csv
				ls -F | grep "/$"|cut -f1 -d'/' > id
				ls -1 */outs/atac_fragments.tsv.gz > frag
				ls -1 */outs/per_barcode_metrics.csv > per
				ls -1 */outs/gex_molecule_info.h5 > gex
				paste id frag per gex -d',' >>aggr.csv
				csv='aggr.csv'
				rm id frag per gex
			else
				echo "Only support atac or multiome data aggregation!"
				exit 1
			fi
		fi
		$cellranger_path aggr --id aggr --csv $csv $ref_type $ref_path --normalize $norm --localcores $threads --localmem $mem
	else
		for i in $(ls -1 ${1}/*R1*.gz);do
			filename=${i##*/}
			prefix=${filename%_S*_L*.gz}
			$cellranger_path count --id $prefix --fastqs $1 $ref_type $ref_path --sample $prefix --localcores $threads --localmem $mem $added_par
		done
	fi
}

main $1

# check running status
if [ $? -ne 0 ]; then
	help
	exit 1
else
	echo "Run succeed"
fi

################ END ################
#          Created by Aone          #
#     quanyiz@stanford.edu      #
################ END ################