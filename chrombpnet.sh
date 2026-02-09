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

# ================= Configuration Summary =================
echo "================================================"
echo "Running ChromBPNet Pipeline"
echo "================================================"
echo "Input List:      $INPUT_FILE"
echo "Genome Fasta:    $GENOME_FA"
echo "Chrom Sizes:     $CHROM_SIZES"
echo "Blacklist:       $BLACKLIST"
echo "Bias Fragments:  ${BIAS_FRAGMENTS:-[Auto Merge]}"
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
echo ">>> Step 1: Running MACS2 Callpeak"

while read -r cluster frag_file; do
    # Skip empty lines or comments
    [[ -z "$cluster" || "$cluster" =~ ^# ]] && continue
    
    echo "  > Calling peaks for ${cluster}..."
    
    if [ ! -f "$frag_file" ]; then
        echo "    [Error] Fragment file not found: $frag_file"
        continue
    fi

    # Check if output exists
    if [ -f "${MACS2_DIR}/${cluster}_peaks.narrowPeak" ]; then
        echo "    [Info] Output exists, skipping MACS2."
        continue
    fi

    # Run in background
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

wait # Wait for all background jobs
echo "Step 1 Done."

# ================= Step 2: Preparing Non-Peaks =================
echo ">>> Step 2: Preparing Non-Peaks (Negatives)"

# 2.1 Slop Blacklist
bedtools slop -i "$BLACKLIST" -g "$CHROM_SIZES" -b 1057 > "${TMP_DIR}/blacklist_slop.bed"

# 2.2 Merge Peaks
echo "  > Merging all peaks..."
# Ensure at least one peak file exists
count_peaks=$(ls "${MACS2_DIR}"/*_peaks.narrowPeak 2>/dev/null | wc -l)
if [ "$count_peaks" -eq 0 ]; then
    echo "Error: No peak files found in ${MACS2_DIR}. Step 1 might have failed."
    exit 1
fi

cat "${MACS2_DIR}"/*_peaks.narrowPeak | grep "^chr" | awk '$8+0 > 2' | sort -k1,1 -k2,2n > "${TMP_DIR}/all_peaks_sorted.bed"

# Copy sorted peaks for Bias training
cp "${TMP_DIR}/all_peaks_sorted.bed" "${MACS2_DIR}/all_peaks.narrowPeak"

# Bedtools Merge & Intersect
bedtools merge -i "${TMP_DIR}/all_peaks_sorted.bed" > "${TMP_DIR}/union_peaks_3col.bed"
intersectBed -a "${TMP_DIR}/union_peaks_3col.bed" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${TMP_DIR}/union_final_3col.bed"

# Format to narrowPeak
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
echo ">>> Step 3: Training Bias Model"

# [Auto-Merge Logic]
if [[ -z "$BIAS_FRAGMENTS" ]]; then
    echo "  [Info] No bias fragments file provided (-b). Automatically merging from input list..."
    
    MERGED_FRAG="${TMP_DIR}/merged_all_fragments.tsv.gz"
    
    # Extract paths from input file
    awk '{print $2}' "$INPUT_FILE" > "${TMP_DIR}/frag_list.txt"
    
    # ------------------ WARNING ------------------
    echo ""
    echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo "  WARNING: Merging large fragment files requires significant"
    echo "           disk space for sorting."
    echo ""
    echo "  Temporary directory: ${TMP_DIR}"
    echo "  Ensure this partition has enough free space (approx. 2-3x total size)."
    echo "  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
    echo ""
    # ---------------------------------------------

    echo "  > Merging, Sorting, and Compressing fragments (This may take a while)..."
    
    # Merge -> Sort (using TMP_DIR) -> Compress
    cat "${TMP_DIR}/frag_list.txt" | xargs zcat | \
    sort -T "$TMP_DIR" -k1,1V -k2,2n | \
    bgzip > "$MERGED_FRAG"
    
    echo "  > Indexing merged fragments..."
    tabix -p bed "$MERGED_FRAG"
    
    BIAS_FRAGMENTS="$MERGED_FRAG"
    echo "  > Merged file created at: $BIAS_FRAGMENTS"
fi

# Validation
if [ ! -f "$BIAS_FRAGMENTS" ]; then
    echo "  [Error] Bias fragments file not found: $BIAS_FRAGMENTS"
    exit 1
fi

# Check if Bias model exists
BIAS_MODEL_PATH=""
if [ -f "${BIAS_DIR}/global/models/bias.h5" ]; then
    BIAS_MODEL_PATH="${BIAS_DIR}/global/models/bias.h5"
elif [ -f "${BIAS_DIR}/global_bias_models_bias.h5" ]; then
    BIAS_MODEL_PATH="${BIAS_DIR}/global_bias_models_bias.h5"
fi

if [ -z "$BIAS_MODEL_PATH" ]; then
    echo "  > Training global bias model using: $(basename "$BIAS_FRAGMENTS")"
    chrombpnet bias pipeline \
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
echo ">>> Step 4: Training Cluster Models"

while read -r cluster frag_file; do
    [[ -z "$cluster" || "$cluster" =~ ^# ]] && continue
    
    echo "Processing Cluster: ${cluster}"
    peak_file="${MACS2_DIR}/${cluster}_peaks.narrowPeak"
    
    for fold in {0..4}; do
        current_out_dir="${MODEL_DIR}/${cluster}/fold_${fold}"
        current_tmp_dir="${TMP_DIR}/${cluster}_fold${fold}"
        mkdir -p "$current_out_dir" "$current_tmp_dir"
        
        # Resume capability
        if [ -f "${current_out_dir}/models/chrombpnet_nobias.h5" ]; then
            echo "    Fold ${fold} already finished. Skipping."
            continue
        fi

        echo "  > Training Fold ${fold}..."
        
        chrombpnet pipeline \
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
            
        echo "    Fold ${fold} Done."
    done
done < "$INPUT_FILE"

echo "=================================================="
echo "All Pipeline Steps Finished Successfully!"
echo "=================================================="