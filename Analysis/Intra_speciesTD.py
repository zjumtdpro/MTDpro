#!/bin/bash
#SBATCH --job-name=TD_analysis
#SBATCH --partition=cpu
#SBATCH -N 1
#SBATCH --ntasks-per-node=1
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --time=336:00:00
#SBATCH --mem-per-cpu=70gb
#SBATCH --cpus-per-task=1

# --------------------------
# Configuration Section
# --------------------------
set -euo pipefail  

REF_ID="" ##Please enter the file address here
SPECIES_NAME="" ##Please enter the file address here
BASE_DIR="" ##Please enter the file address here
REF_FASTA="ref.fna" ##Please enter the file address here

# Path Definitions
DATA_DIR="${BASE_DIR}" ##Please enter the file address here
MUGSY_DIR="${DATA_DIR}/${SPECIES_NAME}" ##Please enter the file address here
OUTPUT_DIR="${MUGSY_DIR}/output" ##Please enter the file address here
TRF_WORKDIR="${BASE_DIR}/" ##Please enter the file address here

# Create directories if missing
mkdir -p "${MUGSY_DIR}" "${OUTPUT_DIR}" "${TRF_WORKDIR}"

# --------------------------
# Main Analysis Pipeline
# --------------------------

# Phase 1: Generate Strain Filter List
filter_list="${MUGSY_DIR}/" ##Please enter the file address here
assembly_list="${DATA_DIR}/rsgb_assembly_summary/" ##Please enter the file address here

> "${filter_list}"  # Initialize output file

while IFS= read -r assembly_id; do
    if [[ "${assembly_id}" == "${REF_ID}" ]]; then
        continue
    fi
    
    genomic_file=$(find "${DATA_DIR}/fungi/" -name "${assembly_id}*genomic.fna.gz" -print -quit)
    
    if [[ -f "${genomic_file}" ]]; then
        echo "${assembly_id}" >> "${filter_list}"
        echo "[INFO] Found assembly: ${assembly_id}"
    else
        echo "[WARNING] Missing assembly: ${assembly_id}" >&2
    fi
done < "${assembly_list}"

# Phase 2: Core Analysis
cd "${MUGSY_DIR}"

process_assembly() {
    local assembly_id=$1
    local assembly_prefix=$(echo "${assembly_id}" | tr '.' '_')
    
    echo "[INFO] Processing assembly: ${assembly_id}"
    
    # File Handling
    cp "${DATA_DIR}/fungi/${assembly_id}*genomic.fna.gz" .
    gzip -d "${assembly_id}.fna.gz"
    
    # MUGSY Alignment
    mugsy --directory "${MUGSY_DIR}" \
          --prefix "${assembly_prefix}" \
          "${REF_FASTA}" \
          "${assembly_id}.fna"
    
    # Process MAF Output
    python "${BASE_DIR}/identify_insertion.py" \
          "${assembly_prefix}.maf" \
          "${REF_ID}" \
          "${assembly_id}_insert.tsv"
    
    # Cleanup
    rm "${assembly_id}.fna" "${assembly_prefix}".*
}

# Parallel processing control
MAX_JOBS=4  # 根据集群资源调整
current_jobs=0

while IFS= read -r assembly_id; do
    process_assembly "${assembly_id}" &
    
    ((current_jobs++))
    if [[ current_jobs -ge MAX_JOBS ]]; then
        wait
        current_jobs=0
    fi
done < "${filter_list}"
wait

# Phase 3: Data Consolidation
consolidate_inserts() {
    # Merge and deduplicate
    cat *_insert.tsv | awk '!seen[$1$2$4]++' > "${SPECIES_NAME}_insert_rmdup.tsv"
    
    # Sequence Extraction
    extract_sequences() {
        local offset=$1
        awk -F'\t' -v off="${offset}" 'BEGIN {OFS="\t"} {
            print $1, $2+off, $2+off+$5
        }' "${SPECIES_NAME}_insert_rmdup.tsv" > tmp_coordinates.bed
        
        bedtools getfasta -fi "${REF_FASTA}" \
                          -bed tmp_coordinates.bed \
                          -fo seq.fasta
        
        awk '/^>/ {getline; print}' seq.fasta | paste "${SPECIES_NAME}_insert_rmdup.tsv" - > combined.tsv
    }
    
    extract_sequences 0
    extract_sequences $5
}

# ... (后续处理保持相似结构，添加错误检查)

# Phase 4: TRF Filtering
run_trf_filter() {
    cd "${TRF_WORKDIR}"
    
    while read -r entry; do
        # ... (添加详细的参数验证)
        trf "${seq_file}" 2 7 7 80 10 30 1000 -f -d -m
        
        python "${BASE_DIR}/trfinder-TRfilter/convert_script2.py"  # 参数化处理
        
        # ... (后续过滤逻辑)
    done < "${OUTPUT_DIR}/before_trf_file.tsv"
}

# Final Statistics
generate_statistics() {
    # ... (数据验证步骤)
    awk -F'\t' 'BEGIN {
        print "Species","AssemblyCount","TotalTDs","UniqueTDs","DivergentStart","ConservedStart2","ConservedStart3"
    } {
        # 统计逻辑...
    }' > "${OUTPUT_DIR}/final_statistics.tsv"
}
