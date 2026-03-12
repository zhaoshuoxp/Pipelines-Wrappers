#!/bin/bash

# ================= Default Configuration =================
# Genome Paths
DEFAULT_GENOME_FA="/home/quanyiz/genome/hg38/hg38.fa"
DEFAULT_CHROM_SIZES="/home/quanyiz/genome/hg38/hg38.chrom.sizes"
DEFAULT_BLACKLIST="/home/quanyiz/genome/hg38/hg38-blacklist.v2.bed"

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

# Downsample Configurations
ENABLE_DOWNSAMPLE=false
TARGET_PEAKS=300000

# ================= Helper Functions =================
usage() {
    echo "Usage: $0 -i <input_list.txt> [options]"
    echo ""
    echo "Required:"
    echo "  -i  Input file (2 columns: ClusterName FragmentPath)"
    echo ""
    echo "Options:"
    echo "  -s  Enable advanced downsampling (Robust Scaled Knee + Proportional Fallback)"
    echo "  -t  Target peak count for downsampling [Default: $TARGET_PEAKS]"
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
while getopts "i:g:c:l:f:d:b:st:h" opt; do
    case "$opt" in
        i) INPUT_FILE="$OPTARG" ;;
        g) GENOME_FA="$OPTARG" ;;
        c) CHROM_SIZES="$OPTARG" ;;
        l) BLACKLIST="$OPTARG" ;;
        f) FOLDS_JSON="$OPTARG" ;;
        d) FOLDS_DIR="$OPTARG" ;;
        b) BIAS_FRAGMENTS="$OPTARG" ;;
        s) ENABLE_DOWNSAMPLE=true ;;
        t) TARGET_PEAKS="$OPTARG" ;;
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

# ================= GPU Detection & Management =================
echo "================================================"
echo ">>> Detecting Available GPUs"
if ! command -v nvidia-smi &> /dev/null; then
    echo "Error: nvidia-smi not found. Cannot detect GPUs."
    exit 1
fi

FREE_MEM_THRESHOLD=2000 

get_hardware_free_gpus() {
    nvidia-smi --query-gpu=index,memory.used --format=csv,noheader,nounits | \
    awk -v thresh="$FREE_MEM_THRESHOLD" -F', ' '$2 < thresh {print $1}'
}

mapfile -t ALL_GPUS < <(nvidia-smi --query-gpu=index --format=csv,noheader)
NUM_GPUS=${#ALL_GPUS[@]}

if [ "$NUM_GPUS" -eq 0 ]; then
    echo "Error: No GPUs available on this system."
    exit 1
fi
echo "Detected $NUM_GPUS total GPU(s) on the system."


# ================= Configuration Summary =================
echo "================================================"
echo "Running ChromBPNet Pipeline"
echo "================================================"
echo "Input List:      $INPUT_FILE"
echo "Genome Fasta:    $GENOME_FA"
echo "Chrom Sizes:     $CHROM_SIZES"
echo "Blacklist:       $BLACKLIST"
echo "Bias Fragments:  ${BIAS_FRAGMENTS:-[Auto Merge]}"
echo "Downsample:      $ENABLE_DOWNSAMPLE (Target: $TARGET_PEAKS)"
echo "Available GPUs:  $NUM_GPUS"
echo "================================================"

# Setup Directory Structure
BASE_DIR=$(pwd)
OUT_DIR="${BASE_DIR}"
MACS2_DIR="${OUT_DIR}/macs2"
BIAS_DIR="${OUT_DIR}/bias"
MODEL_DIR="${OUT_DIR}/models"
TMP_DIR="${OUT_DIR}/tmp"

mkdir -p "$MACS2_DIR" "$BIAS_DIR" "$MODEL_DIR" "$TMP_DIR"
export TMPDIR="${OUT_DIR}/tmp"
export TEMP="${OUT_DIR}/tmp"
export TMP="${OUT_DIR}/tmp"

# ================= Step 1: MACS2 Callpeak =================
echo ">>> Step 1: Running MACS2 Callpeak (CPU Mode)"

while read -r cluster frag_file || [[ -n "$cluster" ]]; do
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

# 2.2 Merge Peaks (Advanced Strategy)
echo "  > Merging all peaks..."
count_peaks=$(ls "${MACS2_DIR}"/*_peaks.narrowPeak 2>/dev/null | wc -l)
if [ "$count_peaks" -eq 0 ]; then
    echo "Error: No peak files found in ${MACS2_DIR}. Step 1 might have failed."
    exit 1
fi

cat "${MACS2_DIR}"/*_peaks.narrowPeak | grep "^chr" | awk '($8+0 > 2) && ($7+0 > 2)' | sort -k1,1 -k2,2n > "${TMP_DIR}/all_peaks_sorted.narrowPeak"
cp "${TMP_DIR}/all_peaks_sorted.narrowPeak" "${MACS2_DIR}/all_peaks.narrowPeak"

if [ "$ENABLE_DOWNSAMPLE" = true ]; then
    echo "  > Advanced Downsampling Enabled. Target: $TARGET_PEAKS"
    
    bedtools cluster -i "${TMP_DIR}/all_peaks_sorted.narrowPeak" | sort -k11,11n -k8,8nr | awk '!seen[$11]++' | cut -f1-10 > "${TMP_DIR}/union_peaks_best.narrowPeak"
    intersectBed -a "${TMP_DIR}/union_peaks_best.narrowPeak" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${MACS2_DIR}/union_baseline.narrowPeak"
    UNION_BASE=$(wc -l < "${MACS2_DIR}/union_baseline.narrowPeak")
    
    if [ "$UNION_BASE" -le "$TARGET_PEAKS" ]; then
        echo "  > Baseline union peaks ($UNION_BASE) is already <= $TARGET_PEAKS. Skipping downsampling."
        cp "${MACS2_DIR}/union_baseline.narrowPeak" "${MACS2_DIR}/union_final.narrowPeak"
    else
        echo "  > Baseline ($UNION_BASE) exceeds $TARGET_PEAKS. Calculating Robust Scaled Composite Score..."
        rm -f "${TMP_DIR}"/*.scored_peaks
        rm -f "${TMP_DIR}/all_knee_peaks.narrowPeak"

        for peak_file in $(ls -1 "${MACS2_DIR}"/*_peaks.narrowPeak|grep -v all_peaks.narrowPeak) ; do
            bname=$(basename "$peak_file")
            awk '($8+0 > 2) && ($7+0 > 2)' "$peak_file" > "${TMP_DIR}/$bname.filtered"
            N=$(wc -l < "${TMP_DIR}/$bname.filtered")
            if [ "$N" -eq 0 ]; then continue; fi

            FC_99=$(awk '{print $7}' "${TMP_DIR}/$bname.filtered" | sort -n | awk '{a[NR]=$1} END {idx=int(NR*0.99); if(idx==0) idx=1; print a[idx]}')
            PVAL_99=$(awk '{print $8}' "${TMP_DIR}/$bname.filtered" | sort -n | awk '{a[NR]=$1} END {idx=int(NR*0.99); if(idx==0) idx=1; print a[idx]}')
            FC_MIN=$(awk '{print $7}' "${TMP_DIR}/$bname.filtered" | sort -n | head -n 1)
            PVAL_MIN=$(awk '{print $8}' "${TMP_DIR}/$bname.filtered" | sort -n | head -n 1)

            awk -v fc99="$FC_99" -v pval99="$PVAL_99" -v fcmin="$FC_MIN" -v pvalmin="$PVAL_MIN" '
            {
                fc = $7; pval = $8;
                if (fc > fc99) fc = fc99;
                if (pval > pval99) pval = pval99;
                
                fcrange = fc99 - fcmin; if(fcrange<=0) fcrange=1e-6;
                pvalrange = pval99 - pvalmin; if(pvalrange<=0) pvalrange=1e-6;
                
                score = ((fc - fcmin) / fcrange) * ((pval - pvalmin) / pvalrange);
                print $0 "\t" score;
            }' "${TMP_DIR}/$bname.filtered" | sort -k11,11nr > "${TMP_DIR}/$bname.scored_peaks"

            KNEE=$(awk '{scores[NR]=$11} END {
                N=NR;
                if(N<3000) { print N; exit; }
                y1=scores[1]; yN=scores[N];
                A = y1 - yN; B = N - 1; C = yN - N * y1;
                max_dist = -1; knee = N;
                for(i=1; i<=N; i++) {
                    dist = A * i + B * scores[i] + C;
                    if(dist < 0) dist = -dist;
                    if(dist > max_dist) { max_dist = dist; knee = i; }
                }
                minkeep = int(N * 0.2);
                if(minkeep < 5000) minkeep = 5000;
                if(minkeep > N) minkeep = N;
                if(knee < minkeep) knee = minkeep;
                print knee;
            }' "${TMP_DIR}/$bname.scored_peaks")

            echo "    - $bname | Total: $N | Knee cut: $KNEE"
            head -n "$KNEE" "${TMP_DIR}/$bname.scored_peaks" >> "${TMP_DIR}/all_knee_peaks.narrowPeak"
        done

        grep "^chr" "${TMP_DIR}/all_knee_peaks.narrowPeak" | sort -k1,1 -k2,2n > "${TMP_DIR}/all_knee_sorted.narrowPeak"

        bedtools cluster -i "${TMP_DIR}/all_knee_sorted.narrowPeak" | \
            sort -k12,12n -k11,11nr | \
            awk '!seen[$12]++' | \
            cut -f1-10 > "${TMP_DIR}/union_knee_best.narrowPeak"

        intersectBed -a "${TMP_DIR}/union_knee_best.narrowPeak" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${MACS2_DIR}/union_knee_final.narrowPeak"

        KNEE_UNION=$(wc -l < "${MACS2_DIR}/union_knee_final.narrowPeak")
        
        if [ "$KNEE_UNION" -ge "$TARGET_PEAKS" ]; then
            echo "  > Knee method yielded >= $TARGET_PEAKS peaks ($KNEE_UNION). Using Knee cutoffs directly!"
            cp "${MACS2_DIR}/union_knee_final.narrowPeak" "${MACS2_DIR}/union_final.narrowPeak"
        else
            echo "  > Knee method yielded $KNEE_UNION < $TARGET_PEAKS. Falling back to Proportional Scaled Score Quota..."
            RATIO=$(awk -v t="$TARGET_PEAKS" -v b="$UNION_BASE" 'BEGIN {print t/b}')
            rm -f "${TMP_DIR}/all_quota_peaks.narrowPeak"
            for scored_file in "${TMP_DIR}"/*.scored_peaks; do
                N=$(wc -l < "$scored_file")
                QUOTA=$(awk -v n="$N" -v r="$RATIO" 'BEGIN {printf "%.0f", n * r}')
                head -n "$QUOTA" "$scored_file" >> "${TMP_DIR}/all_quota_peaks.narrowPeak"
            done

            grep "^chr" "${TMP_DIR}/all_quota_peaks.narrowPeak" | sort -k1,1 -k2,2n > "${TMP_DIR}/all_quota_sorted.narrowPeak"
            
            bedtools cluster -i "${TMP_DIR}/all_quota_sorted.narrowPeak" | \
                sort -k12,12n -k11,11nr | \
                awk '!seen[$12]++' | \
                cut -f1-10 > "${TMP_DIR}/union_quota_best.narrowPeak"
                
            intersectBed -a "${TMP_DIR}/union_quota_best.narrowPeak" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${MACS2_DIR}/union_final.narrowPeak"
            FINAL=$(wc -l < "${MACS2_DIR}/union_final.narrowPeak")
            echo "  > Proportional fallback successful! Final union peaks: $FINAL"
        fi
    fi
else
    echo "  > Downsampling Disabled. Using baseline Cluster & Pick Max strategy."
    bedtools cluster -i "${TMP_DIR}/all_peaks_sorted.narrowPeak" | \
        sort -k11,11n -k8,8nr | \
        awk '!seen[$11]++' | \
        cut -f1-10 > "${TMP_DIR}/union_peaks_best.narrowPeak"
        
    intersectBed -a "${TMP_DIR}/union_peaks_best.narrowPeak" -b "${TMP_DIR}/blacklist_slop.bed" -v > "${MACS2_DIR}/union_final.narrowPeak"
    FINAL=$(wc -l < "${MACS2_DIR}/union_final.narrowPeak")
    echo "  > Final union peaks: $FINAL"
fi

# 2.3 ChromBPNet Prep
if [ ! -d "${MACS2_DIR}/union_filter" ]; then
    echo "  > Running chrombpnet prep nonpeaks..."
    chrombpnet prep nonpeaks \
        -g "$GENOME_FA" \
        -c "$CHROM_SIZES" \
        -p "${MACS2_DIR}/union_final.narrowPeak" \
        -o "${MACS2_DIR}/union_filter" \
        -fl "$FOLDS_JSON" \
        -br "$BLACKLIST" 
else
    echo "  > Nonpeaks already prepared. Skipping."
fi
echo "Step 2 Done."

# ================= Step 3: Training Bias Model =================
echo ">>> Step 3: Training Bias Model (Dynamic GPU Allocation)"

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
    BIAS_GPU=""
    while [ -z "$BIAS_GPU" ]; do
        free_gpus=($(get_hardware_free_gpus))
        if [ ${#free_gpus[@]} -gt 0 ]; then
            BIAS_GPU="${free_gpus[0]}"
        else
            echo "  > [Waiting] No free GPU found for Bias training. Retrying in 30s..."
            sleep 30
        fi
    done

    echo "  > Training global bias model using actual free GPU: ${BIAS_GPU}"
    
    CUDA_VISIBLE_DEVICES=$BIAS_GPU chrombpnet bias pipeline \
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
for gpu in "${ALL_GPUS[@]}"; do
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
            hardware_free_gpus=($(get_hardware_free_gpus))
        
            for gpu in "${hardware_free_gpus[@]}"; do
                pid=${gpu_pids[$gpu]}
                if [ -z "$pid" ] || ! kill -0 "$pid" 2>/dev/null; then
                    assigned_gpu=$gpu
                    break
                fi
            done
            
            if [ -z "$assigned_gpu" ]; then
                wait -n 2>/dev/null || sleep 10
            fi
        done
        
        echo "  > Assigned [${cluster} - Fold ${fold}] to currently available GPU ${assigned_gpu}"
        
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