#!/bin/bash

# ================= Default Configuration =================
# Genome Paths
DEFAULT_GENOME_FA="/home/quanyiz/genome/hg38/hg38.fa"
DEFAULT_CHROM_SIZES="/home/quanyiz/genome/hg38/hg38.chrom.sizes"
DEFAULT_BLACKLIST="/home/quanyiz/genome/hg38/hg38.blacklist.bed"

# Folds Paths
DEFAULT_FOLDS_JSON="/home/quanyiz/genome/bias_models/folds/fold_0.json"
DEFAULT_FOLDS_DIR="/home/quanyiz/genome/bias_models/folds"

# Initialize Variables
INPUT_FILE=""
GENOME_FA="$DEFAULT_GENOME_FA"
CHROM_SIZES="$DEFAULT_CHROM_SIZES"
BLACKLIST="$DEFAULT_BLACKLIST"
FOLDS_JSON="$DEFAULT_FOLDS_JSON"
FOLDS_DIR="$DEFAULT_FOLDS_DIR"
BIAS_FRAGMENTS="" # Defaults to empty, triggering auto-merge logic

# ================= Helper Functions =================
usage() {
    echo "Usage: $0 -i <input_list.txt> [options]"
    echo ""
    echo "Required:"
    echo "  -i  Input file (2 columns: ClusterName FragmentPath)"
    echo ""
    echo "Options:"
    echo "  -b  Bias Fragments File [Default: Auto-merge from input list]"
    echo "  -g  Genome Fasta        [Default: $DEFAULT_GENOME_FA]"
    echo "  -c  Chrom Sizes         [Default: $DEFAULT_CHROM_SIZES]"
    echo "  -l  Blacklist Bed       [Default: $DEFAULT_BLACKLIST]"
    echo "  -f  Fold 0 JSON         [Default: $DEFAULT_FOLDS_JSON]"
    echo "  -d  Folds Directory     [Default: $DEFAULT_FOLDS_DIR]"
    echo "  -h  Show this help message"
    exit 1
}

# ================= Argument Parsing =================
while getopts "i:g:c:l:f:d:b:h" opt; do
    case "$opt" in
        i) INPUT_FILE="$OPTARG" ;;
        g) GENOME_FA="$OPTARG" ;;
        c) CHROM_SIZES="$OPTARG" ;;
        l) BLACKLIST="$OPTARG" ;;
        f) FOLDS_JSON="$OPTARG" ;;
        d) FOLDS_DIR="$OPTARG" ;;
        b) BIAS_FRAGMENTS="$OPTARG" ;;
        h) usage ;;
        *) usage ;;
    esac
done

# Validation
if [[ -z "$INPUT_FILE" ]]; then
    echo "Error: Input file (-i) is required."
    usage
fi

if [ ! -f "$INPUT_FILE" ]; then
    echo "Error: Input file '$INPUT_FILE' not found."
    exit 1
fi

# ================= GPU Detection =================
echo "================================================"
echo ">>> Detecting Available GPUs"
if command -v nvidia-smi &> /dev/null; then
    mapfile -t GPUS < <(nvidia-smi --query-gpu=index --format=csv,noheader)
else
    echo "Error: nvidia-smi not found. Cannot detect GPUs."
    exit 1
fi

NUM_GPUS=${#GPUS[@]}
echo "Detected $NUM_GPUS GPU(s): ${GPUS[*]}"

if [ "$NUM_GPUS" -eq 0 ]; then
    echo "Error: No GPUs available. This pipeline requires GPUs for training."
    exit 1
fi

# ================= Configuration Summary =================
echo "================================================"
echo "Running ChromBPNet Pipeline"
echo "================================================"
echo "Input List:      $INPUT_FILE"
echo "Genome Fasta:    $GENOME_FA"
echo "Chrom Sizes:     $CHROM_SIZES"
echo "Blacklist:       $BLACKLIST"
echo "Bias Fragments:  ${BIAS_FRAGMENTS:-[Auto Merge]}"
echo "Available GPUs:  $NUM_GPUS"
echo "================================================"

# Setup Directory Structure
BASE_DIR=$(pwd)
OUT_DIR="${BASE_DIR}/results"
MACS2_DIR="${OUT_DIR}/macs2"
BIAS_DIR="${OUT_DIR}/bias"
MODEL_DIR="${OUT_DIR}/models"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "$MACS2_DIR" "$BIAS_DIR" "$MODEL_DIR" "$TMP_DIR"

# ================= Step 1: MACS2 Callpeak =================
echo ">>> Step 1: Running MACS2 Callpeak (CPU Mode)"

while read -r cluster frag_file; do
    [[ -z "$cluster" || "$cluster" =~ ^# ]] && continue
    
    echo "  > Calling peaks for ${cluster}..."
    
    if [ ! -f "$frag_file" ]; then
        echo "    [Error] Fragment file not found: $frag_file"
        continue
    fi

    if [ -f "${MACS2_DIR}/${cluster}_peaks.narrowPeak" ]; then
        echo "    [Info] Output exists, skipping MACS2."
        continue
    fi

    macs2 callpeak \
        -t "$frag_file" \
        -f BED \
        -g hs \
        --keep-dup all \
        --pvalue 0.01 \
        --shift -75 \
        --extsize 150 \
        --nomodel \
        --outdir "$MACS2_DIR" \
        -n "${cluster}" > "${MACS2_DIR}/${cluster}.macs2.log" 2>&1 &
        
done < "$INPUT_FILE"

wait # Wait for all background MACS2 jobs
echo "Step 1 Done."

# ================= Step 2: Preparing Non-Peaks =================
echo ">>> Step 2: Preparing Non-Peaks (Negatives)"
# 2.1 Slop Blacklist
bedtools slop -i "$BLACKLIST" -g "$CHROM_SIZES" -b 1057 > "${TMP_DIR}/blacklist_slop.bed"

# 2.2 Merge Peaks
echo "  > Merging all peaks..."
count_peaks=$(ls "${MACS2_DIR}"/*_peaks.narrowPeak 2>/dev/null | wc -l)
if [ "$count_peaks" -eq 0 ]; then
    echo "Error: No peak files found in ${MACS2_DIR}. Step 1 might have failed."
    exit 1
fi

cat "${MACS2_DIR}"/*_peaks.narrowPeak | grep "^chr" | awk '$8+0 > 2' | sort -k1,1 -k2,2n > "${TMP_DIR}/all_peaks_sorted.bed"
cp "${TMP_DIR}/all_peaks_sorted.bed" "${MACS2_DIR}/all_peaks.narrowPeak"

bedtools merge -i "${TMP_DIR}/all_peaks_sorted.bed" > "${TMP_DIR}/union_peaks_3col.bed"
intersectBed -a "${TMP_DIR}/union_peaks_3col.bed" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${TMP_DIR}/union_final_3col.bed"

awk -v OFS="\t" '{
    width = $3 - $2;
    summit = int(width / 2);
    print $1, $2, $3, ".", "0", ".", "0", "0", "0", summit
}' "${TMP_DIR}/union_final_3col.bed" > "${MACS2_DIR}/union_final.narrowPeak"

# 2.3 ChromBPNet Prep
if [ ! -d "${MACS2_DIR}/union_filter" ]; then
    echo "  > Running chrombpnet prep nonpeaks..."
    chrombpnet prep nonpeaks \
        -g "$GENOME_FA" \
        -c "$CHROM_SIZES" \
        -p "${MACS2_DIR}/union_final.narrowPeak" \
        -o "${MACS2_DIR}/union_filter" \
        -fl "$FOLDS_JSON" \
        -br "$BLACKLIST" \
        --tmpdir "$TMP_DIR"
else
    echo "  > Nonpeaks already prepared. Skipping."
fi
echo "Step 2 Done."

# ================= Step 3: Training Bias Model =================
echo ">>> Step 3: Training Bias Model (Using GPU ${GPUS[0]})"

if [[ -z "$BIAS_FRAGMENTS" ]]; then
    echo "  [Info] No bias fragments file provided (-b). Automatically merging from input list..."
    MERGED_FRAG="${TMP_DIR}/merged_all_fragments.tsv.gz"
    awk '{print $2}' "$INPUT_FILE" > "${TMP_DIR}/frag_list.txt"
    
    echo "  > Merging, Sorting, and Compressing fragments (This may take a while)..."
    cat "${TMP_DIR}/frag_list.txt" | xargs zcat | \
    sort -T "$TMP_DIR" -k1,1V -k2,2n | \
    bgzip > "$MERGED_FRAG"
    
    echo "  > Indexing merged fragments..."
    tabix -p bed "$MERGED_FRAG"
    
    BIAS_FRAGMENTS="$MERGED_FRAG"
    echo "  > Merged file created at: $BIAS_FRAGMENTS"
fi

if [ ! -f "$BIAS_FRAGMENTS" ]; then
    echo "  [Error] Bias fragments file not found: $BIAS_FRAGMENTS"
    exit 1
fi

BIAS_MODEL_PATH=""
if [ -f "${BIAS_DIR}/global/models/bias.h5" ]; then
    BIAS_MODEL_PATH="${BIAS_DIR}/global/models/bias.h5"
elif [ -f "${BIAS_DIR}/global_bias_models_bias.h5" ]; then
    BIAS_MODEL_PATH="${BIAS_DIR}/global_bias_models_bias.h5"
fi

if [ -z "$BIAS_MODEL_PATH" ]; then
    echo "  > Training global bias model using: $(basename "$BIAS_FRAGMENTS")"
    
    CUDA_VISIBLE_DEVICES=${GPUS[0]} chrombpnet bias pipeline \
        -ifrag "$BIAS_FRAGMENTS" \
        -d "ATAC" \
        -g "$GENOME_FA" \
        -c "$CHROM_SIZES" \
        -p "${MACS2_DIR}/all_peaks.narrowPeak" \
        -n "${MACS2_DIR}/union_filter_negatives.bed" \
        -fl "$FOLDS_JSON" \
        -b 0.5 \
        -o "$BIAS_DIR/global" \
        --tmpdir "$TMP_DIR"
    
    if [ -f "${BIAS_DIR}/global/models/bias.h5" ]; then
        BIAS_MODEL_PATH="${BIAS_DIR}/global/models/bias.h5"
    else
        echo "Error: Bias model training failed."
        exit 1
    fi
else
    echo "  > Bias model found at $BIAS_MODEL_PATH. Skipping training."
fi

echo "Step 3 Done. Using Bias Model: $BIAS_MODEL_PATH"

# ================= Step 4: Training Cluster Models =================
echo ">>> Step 4: Training Cluster Models (Parallel on $NUM_GPUS GPUs)"

declare -A gpu_pids
for gpu in "${GPUS[@]}"; do
    gpu_pids[$gpu]=""
done

while read -r cluster frag_file; do
    [[ -z "$cluster" || "$cluster" =~ ^# ]] && continue
    
    echo "Queuing Cluster: ${cluster}"
    peak_file="${MACS2_DIR}/${cluster}_peaks.narrowPeak"
    
    for fold in {0..4}; do
        current_out_dir="${MODEL_DIR}/${cluster}/fold_${fold}"
        current_tmp_dir="${TMP_DIR}/${cluster}_fold${fold}"
        mkdir -p "$current_out_dir" "$current_tmp_dir"
        
        if [ -f "${current_out_dir}/models/chrombpnet_nobias.h5" ]; then
            echo "    Fold ${fold} for ${cluster} already finished. Skipping."
            continue
        fi
        
        assigned_gpu=""
        while [ -z "$assigned_gpu" ]; do
            for gpu in "${GPUS[@]}"; do
                pid=${gpu_pids[$gpu]}
                if [ -z "$pid" ] || ! kill -0 "$pid" 2>/dev/null; then
                    assigned_gpu=$gpu
                    break
                fi
            done
            
            if [ -z "$assigned_gpu" ]; then
                wait -n 2>/dev/null || sleep 5
            fi
        done
        
        echo "  > Assigned [${cluster} - Fold ${fold}] to GPU ${assigned_gpu}"
        
        (
            CUDA_VISIBLE_DEVICES=$assigned_gpu chrombpnet pipeline \
                -ifrag "$frag_file" \
                -d "ATAC" \
                -g "$GENOME_FA" \
                -c "$CHROM_SIZES" \
                -p "$peak_file" \
                -n "${MACS2_DIR}/union_filter_negatives.bed" \
                -fl "${FOLDS_DIR}/fold_${fold}.json" \
                -b "$BIAS_MODEL_PATH" \
                -o "$current_out_dir" \
                --tmpdir "$current_tmp_dir" \
                > "${current_out_dir}/train.log" 2>&1
                
            echo "    [Done] ${cluster} - Fold ${fold} finished on GPU ${assigned_gpu}"
        ) &
        
        gpu_pids[$assigned_gpu]=$!
        
    done
done < "$INPUT_FILE"

echo "  > All folds queued. Waiting for remaining background jobs to finish..."
wait

echo "=================================================="
echo "All Pipeline Steps Finished Successfully!"
echo "=================================================="