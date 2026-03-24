#!/bin/bash

# ==========================================
# Default Configurations
# ==========================================
VCF_DIR="/home/quanyiz/General_Use_Data/1000Genome_hg38"
PED_FILE="/home/quanyiz/General_Use_Data/1000Genome_hg38/20130606_g1k_3202_samples_ped_population.txt"
MAF=0.001
WINDOW=500000
MIN_R2=0.2
OUT_SUMMARY="All_Local_Proxy_SNPs.tsv"
INPUT_FILE=""
POP=""

# ==========================================
# Argument Parsing
# ==========================================
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -i|--input) INPUT_FILE="$2"; shift ;;
        -p|--pop) POP="$2"; shift ;;
        -v|--vcf-dir) VCF_DIR="$2"; shift ;;
        -d|--ped) PED_FILE="$2"; shift ;;
        -m|--maf) MAF="$2"; shift ;;
        -w|--window) WINDOW="$2"; shift ;;
        -r|--r2) MIN_R2="$2"; shift ;;
        -o|--out) OUT_SUMMARY="$2"; shift ;;
        -h|--help) 
            echo "Usage: $0 -i <snp_list.txt> -p <POP_CODE> [options]"
            echo "Options:"
            echo "  -v, --vcf-dir   [Path] Directory containing VCFs (Default: $VCF_DIR)"
            echo "  -d, --ped       [Path] PED population file (Default: $PED_FILE)"
            echo "  -m, --maf       [Float] Minimum Allele Frequency (Default: $MAF)"
            echo "  -w, --window    [Int] Search window +/- in bp (Default: $WINDOW)"
            echo "  -r, --r2        [Float] Minimum R-squared threshold (Default: $MIN_R2)"
            echo "  -o, --out       [Path] Output TSV file (Default: $OUT_SUMMARY)"
            exit 0
            ;;
        *) echo "[Error] Unknown parameter passed: $1"; exit 1 ;;
    esac
    shift
done

# ==========================================
# Validation
# ==========================================
if [[ -z "$INPUT_FILE" || -z "$POP" ]]; then
    echo "[Error] Missing required arguments."
    echo "Usage: $0 -i <snp_list.txt> -p <POPULATION>"
    exit 1
fi

if [[ ! -f "$PED_FILE" ]]; then
    echo "[Error] PED file not found: $PED_FILE"
    exit 1
fi

if [[ ! -f "$INPUT_FILE" ]]; then
    echo "[Error] Input SNP list not found: $INPUT_FILE"
    exit 1
fi

# Initialize summary file
echo -e "CHR_A\tBP_A\tTARGET_SNP\tCHR_B\tBP_B\tPROXY_SNP\tR2\tPOPULATION" > "$OUT_SUMMARY"

# ==========================================
# Extract Population Samples
# ==========================================
SAMPLE_LIST="temp_${POP}_samples.txt"
awk -v pop="$POP" 'NR>1 {if ($7 == pop || $6 == pop) print $2}' "$PED_FILE" > "$SAMPLE_LIST"
SAMPLE_COUNT=$(wc -l < "$SAMPLE_LIST")

if [[ "$SAMPLE_COUNT" -eq 0 ]]; then
    echo "[Error] No samples found for population $POP in $PED_FILE"
    rm -f "$SAMPLE_LIST"
    exit 1
fi

echo "[Info] Extracted $SAMPLE_COUNT samples for population $POP."

# ==========================================
# Process SNPs
# ==========================================
while read -r snp; do
    [[ -z "$snp" ]] && continue
    snp=$(echo "$snp" | xargs)

    CHR=""
    POS=""
    USER_TARGET_ID="$snp"

    # 1. Parse SNP format using Regex
    if [[ "$snp" =~ ^rs[0-9]+$ ]]; then
        echo "[Info] rsID detected ($snp). Fetching hg38 coordinates via Ensembl API..."
        res=$(curl -s "https://rest.ensembl.org/variation/human/${snp}?content-type=application/json" | \
              python3 -c "
import sys, json
try:
    d = json.load(sys.stdin)
    m = d.get('mappings', [{}])[0]
    print(m.get('seq_region_name', ''), m.get('start', ''))
except:
    print('','')
" 2>/dev/null)
        CHR=$(echo "$res" | awk '{print $1}')
        POS=$(echo "$res" | awk '{print $2}')
        sleep 0.1
        
    elif [[ "$snp" =~ ^chr([0-9]+|X|Y|MT):([0-9]+) ]]; then
        CHR="${BASH_REMATCH[1]}"
        POS="${BASH_REMATCH[2]}"
    elif [[ "$snp" =~ ^([0-9]+|X|Y|MT):([0-9]+) ]]; then
        CHR="${BASH_REMATCH[1]}"
        POS="${BASH_REMATCH[2]}"
    else
        echo "[Skip] Format not recognized for: $snp"
        continue
    fi

    if [[ -z "$CHR" || -z "$POS" ]]; then
        echo "[Failed] Could not determine coordinates for: $snp"
        continue
    fi

    CHROM="chr${CHR}"
    
    # 2. Auto-scan for correct VCF file
    VCF="${VCF_DIR}/1kGP_high_coverage_Illumina.${CHROM}.filtered.SNV_INDEL_SV_phased_panel.vcf.gz"
    if [[ ! -f "$VCF" ]]; then
        VCF=$(ls ${VCF_DIR}/*${CHROM}.*.vcf.gz ${VCF_DIR}/*${CHROM}_*.vcf.gz 2>/dev/null | head -n 1)
    fi

    if [[ -z "$VCF" || ! -f "$VCF" ]]; then
        echo "[Skip] VCF file not found for chromosome $CHROM in $VCF_DIR"
        continue
    fi

    # 3. Fetch exact Alleles from 1kGP VCF
    ACTUAL_TARGET_ID=$(bcftools query -r ${CHROM}:${POS}-${POS} -f '%CHROM:%POS:%REF:%ALT\n' ${VCF} 2>/dev/null | head -n 1 | tr -d '\r')
    
    if [[ -z "$ACTUAL_TARGET_ID" ]]; then
        echo "[Failed] $USER_TARGET_ID : Position $CHROM:$POS is completely missing in the 1kGP panel."
        continue
    fi

    # 4. Setup PLINK window and variables
    START=$((POS - WINDOW))
    END=$((POS + WINDOW))
    OUT_PREFIX="temp_ld_${POP}_${CHR}_${POS}"
    LD_WINDOW_KB=$(( (WINDOW * 2) / 1000 ))

    # 5. Execute Extraction and PLINK calculation
    bcftools view -r ${CHROM}:${START}-${END} -S ${SAMPLE_LIST} --force-samples ${VCF} -Ou 2>/dev/null | \
    bcftools annotate --set-id '%CHROM:%POS:%REF:%ALT' -Ov -o "${OUT_PREFIX}.vcf" 2>/dev/null
    
    plink --vcf "${OUT_PREFIX}.vcf" \
          --keep-allele-order \
          --allow-extra-chr \
          --r2 \
          --ld-snp "${ACTUAL_TARGET_ID}" \
          --ld-window-r2 ${MIN_R2} \
          --ld-window 99999 \
          --ld-window-kb ${LD_WINDOW_KB} \
          --maf ${MAF} \
          --out "$OUT_PREFIX" > "${OUT_PREFIX}.log" 2>&1
          
    # 6. Parse and Report Results
    if [[ -f "${OUT_PREFIX}.ld" ]]; then
        LINE_COUNT=$(wc -l < "${OUT_PREFIX}.ld")
        if [[ "$LINE_COUNT" -gt 1 ]]; then
            tail -n +2 "${OUT_PREFIX}.ld" | awk -v pop="$POP" -v orig_id="$USER_TARGET_ID" 'BEGIN {OFS="\t"} {print $1,$2,orig_id,$4,$5,$6,$7,pop}' >> "$OUT_SUMMARY"
            PROXIES=$((LINE_COUNT - 1))
            echo "[Success] $USER_TARGET_ID : Found $PROXIES Proxy SNPs (R2 > ${MIN_R2})."
        else
            echo "[Result] $USER_TARGET_ID : No proxies found with R2 > ${MIN_R2}."
        fi
    else
        ERROR_MSG=$(grep "Error:" "${OUT_PREFIX}.log" | head -n 1)
        echo "[Failed] $USER_TARGET_ID : PLINK Error -> ${ERROR_MSG:-Variant might be monomorphic (MAF < $MAF) in $POP}"
    fi
    
    # 7. Cleanup temp files for this SNP
    rm -f "${OUT_PREFIX}.vcf" "${OUT_PREFIX}.log" "${OUT_PREFIX}.nosex" "${OUT_PREFIX}.ld"

done < "$INPUT_FILE"

rm -f "$SAMPLE_LIST"

echo "=================================================="
echo "[Done] Processing complete. Check $OUT_SUMMARY for results."