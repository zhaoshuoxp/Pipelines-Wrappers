#!/bin/bash

# Check dependencies
requires=("cellranger-atac" "cellranger-arc" "cellranger" "cellranger8" "cellranger9")
for i in "${requires[@]}"; do
    if ! command -v "$i" &>/dev/null; then
        echo "Error: $i not found in PATH" >&2
        exit 1
    fi
done

# Default configuration
mem=200
threads=20
norm='none'
aggr='none'
secondary_flag=1  # 1 = default to --nosecondary unless -s specified

# Default reference paths
mm10_path='/home/quanyiz/genome/refdata-cellranger-arc-mm10-2020-A-2.0.0/'
rna_mm10_path='/home/quanyiz/genome/refdata-gex-mm10-2020-A_tdT/'
hg38_path='/home/quanyiz/genome/refdata-cellranger-arc-GRCh38-2020-A-2.0.0/'

# Help message
help(){
	cat <<-EOF
Usage: cellranger.sh [options] <fastq_directory | output_directory>

Options:
  -g [str]    Genome build <hg38|mm10> (required unless -x is set)
  -m [str]    Data type <rna|atac|multiome> (required)
  -x [str]    Custom reference path (overrides -g)
  -t [int]    Threads (default: 20)
  -r [int]    Memory in GB (default: 200)
  -u [int]    Use CellRanger version <8|7> (for -m rna only)
  -a          Run aggregation (aggr) mode
  -c [str]    Custom CSV file for aggr (optional)
  -n          Normalize in aggr mode (default: none)
  -s          Enable secondary analysis (default: off for aggr and multiome)
  --gex_path  RNA fastq path (for multiome mode)
  --atac_path ATAC fastq path (for multiome mode)
  -h          Show help
EOF
	exit 0
}

# If no arguments provided, show help
if [[ $# -eq 0 ]]; then
	help
	exit 1
fi

# Parse arguments (getopt with long options)
TEMP=$(getopt -o g:m:x:t:r:uac:nhs --long gex_path:,atac_path: -n 'cellranger.sh' -- "$@")
eval set -- "$TEMP"

while true; do
	case "$1" in
		--gex_path) gex_path=$2; shift 2 ;;
		--atac_path) atac_path=$2; shift 2 ;;
		-g) genome=$2; shift 2 ;;
		-m)
			if [[ $2 == "rna" ]]; then
				cellranger_path='cellranger9'
				ref_type='--transcriptome'
				added_par='--create-bam true'
				mod='rna'
			elif [[ $2 == "atac" ]]; then
				cellranger_path='cellranger-atac'
				ref_type='--reference'
				mod='atac'
			elif [[ $2 == "multiome" ]]; then
				cellranger_path='cellranger-arc'
				ref_type='--reference'
				mod='arc'
			else
				echo "Only support: rna, atac, multiome"
				exit 1
			fi
			shift 2 ;;
		-x) ref_path=$2; shift 2 ;;
		-t) threads=$2; shift 2 ;;
		-r) mem=$2; shift 2 ;;
		-u)
			if [[ $2 == "8" ]]; then
				cellranger_path='cellranger8'
			elif [[ $2 == "7" ]]; then
				cellranger_path='cellranger'
				added_par=""
			else
				echo "Only version 7 or 8 supported"
				exit 1
			fi
			ref_type='--transcriptome'	
			shift 2 ;;
		-a) aggr='aggr'; shift ;;
		-c) csv=$2; shift 2 ;;
		-n) norm='depth'; shift ;;
		-s) secondary_flag=0; shift ;;  # Explicitly enable secondary analysis
		-h) help ;;
		--) shift; break ;;
		*) echo "Internal error!"; exit 1 ;;
	esac
done

# Select reference from genome
if [[ -z $ref_path && -n $genome ]]; then
	if [[ $genome == "hg38" ]]; then
		ref_path=$hg38_path
	elif [[ $genome == "mm10" ]]; then
		ref_path=$mm10_path
	else
		echo "Unsupported genome. Use hg38 or mm10 or specify with -x"
		exit 1
	fi
fi

# Generate aggregation CSV with absolute paths
generate_aggr_csv() {
	local mode=$1
	local out="aggr.csv"
	echo "Generating $out for $mode..."

	if [[ $mode == "atac" ]]; then
		echo "library_id,fragments,cells" > "$out"
		find -L . -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
			id=$(basename "$dir")
			frag="$dir/outs/fragments.tsv.gz"
			cell="$dir/outs/singlecell.csv"
			[[ -f $frag && -f $cell ]] && echo "$id,$(realpath "$frag"),$(realpath "$cell")" >> "$out"
		done
	elif [[ $mode == "arc" ]]; then
		echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > "$out"
		find -L . -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
			id=$(basename "$dir")
			frag="$dir/outs/atac_fragments.tsv.gz"
			per="$dir/outs/per_barcode_metrics.csv"
			gex="$dir/outs/gex_molecule_info.h5"
			[[ -f $frag && -f $per && -f $gex ]] && echo "$id,$(realpath "$frag"),$(realpath "$per"),$(realpath "$gex")" >> "$out"
		done
	fi
	csv="$out"
}

# Main execution
main() {
	if [[ -z $mod ]]; then echo "Missing -m (data type)"; exit 1; fi
	if [[ -z $ref_path ]]; then echo "Missing -g or -x (reference)"; exit 1; fi

	# Determine final value of --nosecondary
	secondary=""
	if [[ $secondary_flag -eq 1 && ( $mod == "arc" || $aggr == "aggr" ) ]]; then
		secondary="--nosecondary"
	fi

	# MULTIOME MODE
	if [[ $mod == "arc" && $aggr != "aggr" ]]; then
		if [[ -z $gex_path || -z $atac_path ]]; then
			echo "Error: --gex_path and --atac_path required in multiome mode"
			exit 1
		fi

		echo "🔍 Matching RNA and ATAC fastqs by core sample name..."

		declare -A rna_prefix_map
		declare -A atac_prefix_map

		while read -r f; do
			fname=$(basename "$f")
			core=$(echo "$fname" | sed -E 's/^([^-_]+[-_][0-9]+).*$/\1/')
			prefix="${fname%%_S*}"
			rna_prefix_map["$core"]="$prefix"
		done < <(find "$gex_path" -name "*_R1*.f*q.gz")

		while read -r f; do
			fname=$(basename "$f")
			core=$(echo "$fname" | sed -E 's/^([^-_]+[-_][0-9]+).*$/\1/')
			prefix="${fname%%_S*}"
			atac_prefix_map["$core"]="$prefix"
		done < <(find "$atac_path" -name "*_R1*.f*q.gz")

		abs_gex_path=$(realpath "$gex_path")
		abs_atac_path=$(realpath "$atac_path")

		for core in "${!rna_prefix_map[@]}"; do
			if [[ -n "${atac_prefix_map[$core]}" ]]; then
				rna_prefix=${rna_prefix_map[$core]}
				atac_prefix=${atac_prefix_map[$core]}
				echo "🧬 Running: $core → RNA=$rna_prefix  ATAC=$atac_prefix"

				arc_csv="arc_input_${core}.csv"
				echo "fastqs,sample,library_type" > "$arc_csv"
				echo "$abs_gex_path,$rna_prefix,Gene Expression" >> "$arc_csv"
				echo "$abs_atac_path,$atac_prefix,Chromatin Accessibility" >> "$arc_csv"

				$cellranger_path count \
					--id "${core}_arc" \
					--libraries "$arc_csv" \
					--reference "$ref_path" \
					--localcores "$threads" \
					--localmem "$mem" \
					$secondary
			else
				echo "⚠️  No matching ATAC for RNA sample $core"
			fi
		done
		return
	fi

	# AGGR mode
	if [[ $aggr == "aggr" ]]; then
		[[ -z $csv ]] && generate_aggr_csv "$mod"
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
			echo "🧬 Processing: $prefix"
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

main "$1"

if [[ $? -eq 0 ]]; then
	echo "✅ CellRanger run completed successfully."
else
	echo "❌ CellRanger run failed."
	help
	exit 1
fi
