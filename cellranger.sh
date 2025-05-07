#!/bin/bash

# check dependencies
requires=("cellranger-atac" "cellranger-arc" "cellranger" "cellranger8" "cellranger9")
for i in "${requires[@]}"; do
    if ! command -v "$i" &>/dev/null; then
        echo "Error: $i not found in PATH" >&2
        exit 1
    fi
done

# default settings
mem=200
threads=20
norm='none'
aggr='none'
secondary='--nosecondary'

# default reference paths
mm10_path='/home/quanyiz/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0/'
rna_mm10_path='/home/quanyiz/genome/refdata-gex-mm10-2020-A_tdT/'
hg38_path='/home/quanyiz/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/'

# help message
help(){
	cat <<-EOF
Usage: cellranger.sh [options] <fastq_directory | output_directory>

Options:
  -g [str] Genome build <hg38|mm10> (required if no -x)
  -m [str] Data type <rna|atac|multiome> (required)
  -x [str] Custom reference path (overrides -g)
  -t [int] Threads (default: 20)
  -r [int] Memory in GB (default: 200)
  -u [int] Use specific CellRanger version <8|7> (only for -m rna)
  -a       Run "aggr" mode instead of "count"
  -c [str] Custom CSV for aggr (only with -a)
  -n       Normalize in aggr mode (default: none)
  -s       Enable secondary analysis (default: off)
  -h       Show help
EOF
	exit 0
}

# parse arguments
while getopts "g:m:x:t:r:uac:nsh" arg; do
	case $arg in
		g)
			genome=$OPTARG
			if [[ $genome == "hg38" ]]; then
				ref_path=$hg38_path
			elif [[ $genome == "mm10" ]]; then
				ref_path=$mm10_path
			else
				echo "Only hg38 or mm10 supported unless using -x"
				exit 1
			fi
			;;
		m)
			if [[ $OPTARG == "rna" ]]; then
				cellranger_path='cellranger9'
				ref_type='--transcriptome'
				added_par='--create-bam true'
				mod='rna'
			elif [[ $OPTARG == "atac" ]]; then
				cellranger_path='cellranger-atac'
				ref_type='--reference'
				mod='atac'
			elif [[ $OPTARG == "multiome" ]]; then
				cellranger_path='cellranger-arc'
				ref_type='--reference'
				mod='arc'
			else
				echo "Only support: rna, atac, multiome"
				exit 1
			fi
			;;
		x) ref_path=$OPTARG ;;
		t) threads=$OPTARG ;;
		r) mem=$OPTARG ;;
		u)
			if [[ $OPTARG == "8" ]]; then
				cellranger_path='cellranger8'
			elif [[ $OPTARG == "7" ]]; then
				cellranger_path='cellranger'
				added_par=""
			else
				echo "Only version 7 or 8 supported"
				exit 1
			fi
			ref_type='--transcriptome'
			;;
		a) aggr='aggr' ;;
		c) csv=$OPTARG ;;
		n) norm='depth' ;;
		s) secondary='' ;;
		h) help ;;
		*) help; exit 1 ;;
	esac
done

shift $((OPTIND - 1))

# generate aggr CSV
generate_aggr_csv() {
	local mode=$1
	local out="aggr.csv"
	echo "Generating $out for $mode..."
	if [[ $mode == "atac" ]]; then
		echo "library_id,fragments,cells" > "$out"
		find . -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
			id=$(basename "$dir")
			frag="$dir/outs/fragments.tsv.gz"
			cell="$dir/outs/singlecell.csv"
			[[ -f $frag && -f $cell ]] && echo "$id,$(realpath $frag),$(realpath $cell)" >> "$out"
		done
	elif [[ $mode == "arc" ]]; then
		echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > "$out"
		find . -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
			id=$(basename "$dir")
			frag="$dir/outs/atac_fragments.tsv.gz"
			per="$dir/outs/per_barcode_metrics.csv"
			gex="$dir/outs/gex_molecule_info.h5"
			[[ -f $frag && -f $per && -f $gex ]] && echo "$id,$(realpath $frag),$(realpath $per),$(realpath $gex)" >> "$out"
		done
	fi
	csv="$out"
}

# main function
main() {
	if [[ -z $mod ]]; then echo "Missing -m (data type)"; exit 1; fi
	if [[ -z $ref_path ]]; then echo "Missing -g or -x (reference)"; exit 1; fi

	# aggr mode
	if [[ $aggr == "aggr" ]]; then
		if [[ -z $csv ]]; then
			generate_aggr_csv "$mod"
		fi
		$cellranger_path aggr \
			--id aggr_$(date +%Y%m%d) \
			--csv "$csv" \
			$ref_type "$ref_path" \
			--normalize "$norm" \
			--localcores "$threads" \
			--localmem "$mem" \
			$secondary
	else
		if [[ $mod == "rna" && $genome == "mm10" ]]; then
			ref_path=$rna_mm10_path
		fi
		for i in "$1"/*_R1*.f*q.gz; do
			[[ ! -e $i ]] && continue
			filename=$(basename "$i")
			prefix="${filename%%_S*}"
			echo "Running sample: $prefix"
			$cellranger_path count \
				--id "$prefix" \
				--fastqs "$1" \
				$ref_type "$ref_path" \
				--sample "$prefix" \
				--localcores "$threads" \
				--localmem "$mem" \
				$added_par $secondary
		done
	fi
}

# run main
main "$1"

# report
if [[ $? -eq 0 ]]; then
	echo "✅ CellRanger run completed successfully."
else
	echo "❌ CellRanger run failed."
	help
	exit 1
fi
