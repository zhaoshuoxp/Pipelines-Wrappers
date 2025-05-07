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
secondary='--nosecondary'

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
  -s          Enable secondary analysis (default: off)
  --gex_path  RNA fastq path (required for multiome)
  --atac_path ATAC fastq path (required for multiome)
  -h          Show help
EOF
	exit 0
}

# Parse arguments (getopt with long options)
TEMP=$(getopt -o g:m:x:t:r:uac:nh --long gex_path:,atac_path: -n 'cellranger.sh' -- "$@")
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
		-s) secondary=''; shift ;;
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

# Generate aggregation CSV
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
			[[ -f $frag && -f $cell ]] && echo "$id,$(realpath "$frag"),$(realpath "$cell")" >> "$out"
		done
	elif [[ $mode == "arc" ]]; then
		echo "library_id,atac_fragments,per_barcode_metrics,gex_molecule_info" > "$out"
		find . -mindepth 1 -maxdepth 1 -type d | while read -r dir; do
			id=$(basename "$dir")
			frag="$dir/outs/atac_fragments.tsv.gz"
			per="$dir/outs/per_barcode_metrics.csv"
			gex="$dir/outs/gex_molecule_info.h5"
			[[ -f $frag && -f $per && -f $gex ]] && echo "$id,$(realpath "$frag"),$(realpath "$per"),$(realpath "$gex")" >> "$out"
		done
	fi
	csv="$out"
}

# Main execution logic
main() {
	if [[ -z $mod ]]; then echo "Missing -m (data type)"; exit 1; fi
	if [[ -z $ref_path ]]; then echo "Missing -g or -x (reference)"; exit 1; fi

	# Multiome mode with loose core name matching
	if [[ $mod == "arc" && $aggr != "aggr" ]]; then
		if [[ -z $gex_path || -z $atac_path ]]; then
			echo "Error: --gex_path and --atac_path required in multiome mode"
			exit 1
		fi

		echo "🔍 Matching RNA and ATAC fastqs by core sample name..."

		declare -A rna_map
		declare -A atac_map

		# RNA: extract core like Mult_1 from filename
		while read -r f; do
			fname=$(basename "$f")
			core=$(echo "$fname" | sed -E 's/_RNA.*$//' | sed -E 's/[-_][^_]+_S[0-9]+_L[0-9]+_R[12]_001\.f.*gz//' )
			rna_map["$core"]=$gex_path
		done < <(find "$gex_path" -name "*_R1*.f*q.gz")

		while read -r f; do
			fname=$(basename "$f")
			core=$(echo "$fname" | sed -E 's/[-_][^_]+_S[0-9]+_L[0-9]+_R[12]_001\.f.*gz//' )
			atac_map["$core"]=$atac_path
		done < <(find "$atac_path" -name "*_R1*.f*q.gz")

		for sample in "${!rna_map[@]}"; do
			if [[ -n "${atac_map[$sample]}" ]]; then
				echo "🧬 Running: $sample"
				arc_csv="arc_input_${sample}.csv"
				echo "fastqs,sample,library_type" > "$arc_csv"
				echo "${rna_map[$sample]},$sample,Gene Expression" >> "$arc_csv"
				echo "${atac_map[$sample]},$sample,Chromatin Accessibility" >> "$arc_csv"

				$cellranger_path count \
					--id "${sample}_arc" \
					--libraries "$arc_csv" \
					--reference "$ref_path" \
					--localcores "$threads" \
					--localmem "$mem" \
					$secondary
			else
				echo "⚠️  No matching ATAC for RNA sample $sample"
			fi
		done
		return
	fi

	# Aggregation mode
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

# Run it
main "$1"

if [[ $? -eq 0 ]]; then
	echo "✅ CellRanger run completed successfully."
else
	echo "❌ CellRanger run failed."
	help
	exit 1
fi
