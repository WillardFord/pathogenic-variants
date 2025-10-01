#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DATA_DIR="${PROJECT_ROOT}/data"
OUTPUT_DIR="${PROJECT_ROOT}/output"
TMP_DIR="${PROJECT_ROOT}/tmp"
DOWNLOAD_LIST_PATH="${OUTPUT_DIR}/vcf_download_list.txt"
OUTPUT_DIR="${PROJECT_ROOT}/output/gene_filtered"

CLINVAR_BASE_PATH="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf/"
CLINVAR_ANNOTATION="${DATA_DIR}/clinvar_with_chr.vcf.gz"
CLINVAR_ANNOTATION_TBI="${CLINVAR_ANNOTATION}.tbi"
GENE_LIST="${DATA_DIR}/pathogenic_genes.txt"
GENE_BED="${DATA_DIR}/pathogenic_genes_1000000bp.bed"
GENE_INTERVAL_LIST="${DATA_DIR}/pathogenic_genes_1000000bp.interval_list"
START_INDEX="${START_INDEX:-0}"
BATCH_SIZE=300
GSUTIL_THREADS=32

readonly CLINVAR_BASE_PATH CLINVAR_ANNOTATION CLINVAR_ANNOTATION_TBI GENE_LIST GENE_BED GENE_INTERVAL_LIST START_INDEX BATCH_SIZE GSUTIL_THREADS

for path in "${CLINVAR_ANNOTATION}" "${CLINVAR_ANNOTATION_TBI}" "${GENE_LIST}" "${GENE_BED}" "${GENE_INTERVAL_LIST}"; do
  if [[ ! -f "${path}" ]]; then
    echo "Missing required file: ${path}" >&2
    exit 1
  fi
done

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"
GENE_PATTERN=$(paste -sd'|' "${GENE_LIST}")
if [[ -z "${GENE_PATTERN}" ]]; then
  echo "Gene list produced an empty pattern" >&2
  exit 1
fi
FILTER_EXPR="INFO/GENEINFO ~ \"(${GENE_PATTERN})\""

echo "Using filter expression: ${FILTER_EXPR}"
echo "Using gene interval list: ${GENE_INTERVAL_LIST}"

echo "START_INDEX set to ${START_INDEX}"

# Check if we already have a download list - if so, skip interval processing entirely
if [[ -f "${DOWNLOAD_LIST_PATH}" ]]; then
  echo "Found existing download list at ${DOWNLOAD_LIST_PATH}; skipping interval processing."
else
  echo "Listing interval shards..."
  mapfile -t interval_files < <(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.interval_list 2>/dev/null)

  if (( ${#interval_files[@]} > 0 )); then
    echo "Found ${#interval_files[@]} interval_list file(s); screening in batches of ${BATCH_SIZE}."
    excluded_intervals=0
    declare -A shard_selected=()

    # Start a fresh persistent download list and append as shards are discovered
    : > "${DOWNLOAD_LIST_PATH}"

    total_intervals=${#interval_files[@]}
    if (( START_INDEX >= total_intervals )); then
      echo "START_INDEX (${START_INDEX}) exceeds available interval shards (${total_intervals}); nothing to filter."
      exit 0
    fi

    # Iterate over all interval files to determine which vcfs have overlaps with gene interval list.
    index=${START_INDEX}
    while (( index < total_intervals )); do
      batch_paths=("${interval_files[@]:index:BATCH_SIZE}")
      index=$((index + ${#batch_paths[@]}))

      batch_dir=$(mktemp -d "${TMP_DIR}/interval_batch.XXXXXX")
      batch_list=$(mktemp "${TMP_DIR}/interval_batch_list.XXXXXX.txt")
      printf '%s\n' "${batch_paths[@]}" > "${batch_list}"

      echo "Downloading interval batch to ${batch_dir} (size=${#batch_paths[@]})"
      gsutil -u "${GOOGLE_PROJECT}" -m \
        -o "GSUtil:parallel_thread_count=${GSUTIL_THREADS}" \
        -o "GSUtil:parallel_process_count=1" \
        cp -I "${batch_dir}/" < "${batch_list}"
      rm -f "${batch_list}"

      combined_file=$(mktemp "${TMP_DIR}/batch.XXXXXX.interval_list")
      overlap_file=$(mktemp "${TMP_DIR}/batch_overlap.XXXXXX.interval_list")
      header_written=false

      for interval_path in "${batch_paths[@]}"; do
        prefix=$(basename "${interval_path}" .interval_list)
        local_interval="${batch_dir}/${prefix}.interval_list"
        if [[ ! -f "${local_interval}" ]]; then
          echo "Warning: interval ${local_interval} missing after download; skipping shard ${prefix}" >&2
          excluded_intervals=$((excluded_intervals + 1))
          continue
        fi

        # Prepend the number in name of each individual interval to combined file so we can track which interval is which.
        if [[ ${header_written} == false ]]; then
          grep '^@' "${local_interval}" >> "${combined_file}"
          header_written=true
        fi

        awk -v shard="${prefix}" 'BEGIN{OFS="\t"} !/^@/ { $5 = shard; print }' "${local_interval}" >> "${combined_file}"
      done

      if [[ ${header_written} == false ]]; then
        rm -f "${combined_file}" "${overlap_file}"
        rm -rf "${batch_dir}"
        continue
      fi

      # Intersection?
      gatk IntervalListTools \
        -I "${combined_file}" \
        -SI "${GENE_INTERVAL_LIST}" \
        --ACTION OVERLAPS \
        -O "${overlap_file}" \
        --QUIET true

      # If yes, then add corresponding vcfs to to_process list
      if tail -n +2 "${overlap_file}" | grep -v '^@' | grep -q .; then
        while IFS= read -r shard_id; do
          [[ -z "${shard_id}" ]] && continue
          if [[ -z "${shard_selected[${shard_id}]:-}" ]]; then
            shard_selected["${shard_id}"]=1
            # Append discovered shard files to the persistent list immediately
            printf '%s\n' "${CLINVAR_BASE_PATH}${shard_id}.vcf.bgz" >> "${DOWNLOAD_LIST_PATH}"
            printf '%s\n' "${CLINVAR_BASE_PATH}${shard_id}.vcf.bgz.tbi" >> "${DOWNLOAD_LIST_PATH}"
          fi
        done < <(tail -n +2 "${overlap_file}" | grep -v '^@' | awk -F'\t' '{print $5}' | sort -u)
      else
        excluded_intervals=$((excluded_intervals + ${#batch_paths[@]}))
      fi

      rm -f "${combined_file}" "${overlap_file}"
      rm -rf "${batch_dir}"
    done

    echo "Selected shard(s) after interval screening (skipped ${excluded_intervals})."
  fi
fi

process_start=0

# Download list should always exist at this point
if [[ ! -f "${DOWNLOAD_LIST_PATH}" ]]; then
  echo "Download list not found at ${DOWNLOAD_LIST_PATH}. This should not happen." >&2
  exit 1
fi

# We'll get the actual count from the download list
total=$(awk '{print $0}' "${DOWNLOAD_LIST_PATH}" | sed -E 's#.*/([^/]+)\.vcf\.bgz(\.tbi)?$#\1#' | grep -E '^[0-9]+' | sort -u | wc -l)
last_index=$((total - 1))

processed=${process_start}
processed_with_times=0
skipped=0
copy_time_total=0
filter_time_total=0

format_avg() {
  local total_time=$1
  local count=$2
  if (( count == 0 )); then
    printf 'n/a'
  else
    awk -v t="$total_time" -v c="$count" 'BEGIN { printf "%.2f", t / c }'
  fi
}


# Build a list of shards that still need processing and download them in bulk
bulk_dir=$(mktemp -d "${TMP_DIR}/vcf_bulk.XXXXXX")
download_list_file="${DOWNLOAD_LIST_PATH}"
to_process_prefixes=()

# Always derive prefixes from the persistent download list
if [[ -f "${download_list_file}" ]]; then
  echo "Using download list at ${download_list_file}."
  mapfile -t to_process_prefixes < <(awk '{print $0}' "${download_list_file}" | sed -E 's#.*/([^/]+)\.vcf\.bgz(\.tbi)?$#\1#' | grep -E '^[0-9]+' | sort -u)
else
  echo "Download list not found at ${download_list_file}. Interval screening should have created it." >&2
  rm -rf "${bulk_dir}"
  exit 1
fi

# Check if any outputs are already completed and skip them
completed_prefixes=()
for prefix in "${to_process_prefixes[@]}"; do
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"
  output_tbi="${output_file}.tbi"
  
  if [[ -f "${output_file}" && -f "${output_tbi}" ]]; then
    echo "Skipping shard ${prefix}: already processed."
    skipped=$((skipped + 1))
    completed_prefixes+=("${prefix}")
  fi
done

# Remove completed prefixes from the processing list
for completed in "${completed_prefixes[@]}"; do
  to_process_prefixes=("${to_process_prefixes[@]/$completed}")
done

if (( ${#to_process_prefixes[@]} == 0 )); then
  echo "Nothing new to process; all shard(s) already completed."
  rm -rf "${bulk_dir}"
  avg_copy=$(format_avg "$copy_time_total" "$processed_with_times")
  avg_filter=$(format_avg "$filter_time_total" "$processed_with_times")
  echo "Completed: ${processed} shard(s); skipped ${skipped}."
  echo "Average copy time (processed shards): ${avg_copy}s"
  echo "Average filter time (processed shards): ${avg_filter}s"
  exit 0
fi

echo "Bulk downloading ${#to_process_prefixes[@]} shard(s) to ${bulk_dir} using parallel gsutil..."
copy_start=$(date +%s)
gsutil -u "${GOOGLE_PROJECT}" -m \
  -o "GSUtil:parallel_thread_count=${GSUTIL_THREADS}" \
  -o "GSUtil:parallel_process_count=1" \
  cp -I "${bulk_dir}/" < "${download_list_file}"
copy_end=$(date +%s)
copy_time_total=$((copy_time_total + (copy_end - copy_start)))

# Process downloaded VCFs in batches with parallelization
PROCESS_BATCH_SIZE=30
MAX_PARALLEL_BATCHES=8  # Number of batches to process in parallel
total_prefixes=${#to_process_prefixes[@]}
batch_count=0

# Function to process a single batch
process_batch() {
  local batch_num=$1
  local start_idx=$2
  local end_idx=$3
  local batch_prefixes=("${@:4}")
  
  echo "Processing batch ${batch_num}: ${#batch_prefixes[@]} files (${start_idx}-${end_idx} of ${total_prefixes})"
  
  # Create list of VCF files for this batch
  batch_vcf_list=$(mktemp "${TMP_DIR}/batch_vcfs.${batch_num}.XXXXXX.txt")
  for prefix in "${batch_prefixes[@]}"; do
    local_vcf="${bulk_dir}/${prefix}.vcf.bgz"
    if [[ -f "${local_vcf}" ]]; then
      echo "${local_vcf}" >> "${batch_vcf_list}"
    else
      echo "Warning: missing ${local_vcf}; skipping from batch ${batch_num}." >&2
    fi
  done
  
  # Check if we have any files in this batch
  if [[ ! -s "${batch_vcf_list}" ]]; then
    echo "No valid files in batch ${batch_num}; skipping."
    rm -f "${batch_vcf_list}"
    return
  fi
  
  # Concat the batch (faster than merge for same samples)
  concat_vcf="${OUTPUT_DIR}/batched_vcfs/batch_${batch_num}.vcf.bgz"
  mkdir -p "${OUTPUT_DIR}/batched_vcfs"
  echo "Concatenating ${#batch_prefixes[@]} VCF files for batch ${batch_num}..."
  bcftools concat -f "${batch_vcf_list}" -Oz -o "${concat_vcf}"
  tabix -p vcf "${concat_vcf}"
  rm -f "${batch_vcf_list}"
  
  # Annotate and filter the concatenated file
  bc_start=$(date +%s)
  bcftools annotate \
    -a "${CLINVAR_ANNOTATION}" \
    -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
    -k \
    -Ou "${concat_vcf}" | \
  bcftools view \
    -i "${FILTER_EXPR}" \
    -Oz -o "${OUTPUT_DIR}/batch_${batch_num}_filtered.vcf.bgz"
  tabix -p vcf "${OUTPUT_DIR}/batch_${batch_num}_filtered.vcf.bgz"
  bc_end=$(date +%s)
  
  echo "Filtered batch ${batch_num}: ${#batch_prefixes[@]} files to ${OUTPUT_DIR}/batch_${batch_num}_filtered.vcf.bgz"
  echo "Batched VCF saved to: ${concat_vcf}"
}

# Process batches in parallel
for (( i=0; i<total_prefixes; i+=PROCESS_BATCH_SIZE )); do
  batch_count=$((batch_count + 1))
  batch_end=$((i + PROCESS_BATCH_SIZE))
  if (( batch_end > total_prefixes )); then
    batch_end=${total_prefixes}
  fi
  
  batch_prefixes=("${to_process_prefixes[@]:i:batch_end-i}")
  
  # Start batch processing in background
  process_batch "${batch_count}" "${i}" "${batch_end}" "${batch_prefixes[@]}" &
  
  # Limit number of parallel jobs
  if (( batch_count % MAX_PARALLEL_BATCHES == 0 )); then
    wait  # Wait for all background jobs to complete
  fi
done

# Wait for any remaining background jobs
wait

# Clean up all bulk-downloaded files
rm -rf "${bulk_dir}"

avg_copy=$(format_avg "$copy_time_total" "$processed_with_times")
avg_filter=$(format_avg "$filter_time_total" "$processed_with_times")

echo "Completed: ${processed}/${total} shard(s); skipped ${skipped}."
echo "Average copy time (processed shards): ${avg_copy}s"
echo "Average filter time (processed shards): ${avg_filter}s"
