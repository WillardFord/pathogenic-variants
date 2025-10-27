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

[[ -d "${INPUT_DIR}" ]] || { echo "Error: ${INPUT_DIR} not found"; exit 1; }
[[ -f "${GENE_INTERVAL_LIST}" ]] || { echo "Error: ${GENE_INTERVAL_LIST} not found"; exit 1; }

mkdir -p "${OUTPUT_DIR}"

# Load all interval_list files into an array
mapfile -t interval_files < <(find "${INPUT_DIR}" -type f -name '*.interval_list')

# Process all interval_list files in parallel
export OUTPUT_DIR GENE_INTERVAL_LIST
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
done

# Wait for all background jobs to complete
wait

echo "Completed processing interval lists"