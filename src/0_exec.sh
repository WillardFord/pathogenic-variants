#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)

SCRIPTS=(
  "4.5_flatten_vcfs.sh"
  "5_label_vcfs.sh"
  "6_filter_vcf.sh"
)

for script in "${SCRIPTS[@]}"; do
  echo "Running ${script}..."
  "${SCRIPT_DIR}/${script}"
done

echo "Running vcf_to_parquet conversion..."
uv run python "${SCRIPT_DIR}/../analysis/vcf_to_parquet.py"

echo "All steps completed."
