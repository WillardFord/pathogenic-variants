#!/usr/bin/env bash
set -euo pipefail

# Select overlapping intervals from the gene interval list
# All inputs lie in output/interval_lists/
# If output of `grep "^chr"` is not empty then we should include the file
# Else we can skip it

mkdir -p output
gs_base="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf"
output_file="output/vcf_download_list.txt"
echo -n > "${output_file}"

interval_list_dir="output/interval_lists_with_overlapping_ranges"

# Read file list into array to avoid glob expansion issues
mapfile -t files < <(ls "${interval_list_dir}"/*.interval_list 2>/dev/null)


for file in "${files[@]}"; do
  if grep -m 1 -q "^chr" "$file"; then
    echo "${gs_base}/${file%.interval_list}.vcf.bgz" >> "${output_file}"
  fi
done

echo "Wrote ${#files[@]} interval lists to ${output_file}"