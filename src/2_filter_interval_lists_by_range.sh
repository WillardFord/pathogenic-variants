#!/usr/bin/env bash
set -euo pipefail

# Simple, highly parallel script to find overlap between interval lists and gene intervals
# Input: output/interval_lists/*.interval_list
# Output: output/interval_lists_with_overlapping_ranges/*.interval_list

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

INPUT_DIR="${PROJECT_ROOT}/output/interval_lists"
OUTPUT_DIR="${PROJECT_ROOT}/output/interval_lists_with_overlapping_ranges"
GENE_INTERVAL_LIST="${PROJECT_ROOT}/data/pathogenic_genes_1000000bp.interval_list"
MAX_JOBS="${MAX_JOBS:-$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 8)}"

[[ -d "${INPUT_DIR}" ]] || { echo "Error: ${INPUT_DIR} not found"; exit 1; }
[[ -f "${GENE_INTERVAL_LIST}" ]] || { echo "Error: ${GENE_INTERVAL_LIST} not found"; exit 1; }

mkdir -p "${OUTPUT_DIR}"

# Load all interval_list files into an array
mapfile -t interval_files < <(find "${INPUT_DIR}" -type f -name '*.interval_list')

# Process all interval_list files in parallel (with job limit)
export OUTPUT_DIR GENE_INTERVAL_LIST
active_jobs=0
for interval_file in "${interval_files[@]}"; do
  {
    output_file="${OUTPUT_DIR}/$(basename "${interval_file}")"
    gatk IntervalListTools \
      -I "${interval_file}" \
      -SI "${GENE_INTERVAL_LIST}" \
      --ACTION OVERLAPS \
      -O "${output_file}" \
      --QUIET true
  } &
  
  active_jobs=$((active_jobs + 1))
  
  # Wait when we hit the job limit
  if (( active_jobs >= MAX_JOBS )); then
    wait -n  # Wait for any background job to finish
    active_jobs=$((active_jobs - 1))
  fi
done

# Wait for all remaining background jobs to complete
wait

echo "Completed processing interval lists"