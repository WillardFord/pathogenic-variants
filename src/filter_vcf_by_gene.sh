#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DATA_DIR="${PROJECT_ROOT}/data"
OUTPUT_DIR="${PROJECT_ROOT}/output"
TMP_DIR="${PROJECT_ROOT}/tmp"

CLINVAR_BASE_PATH="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf/"
CLINVAR_ANNOTATION="${DATA_DIR}/clinvar.vcf.gz"
CLINVAR_ANNOTATION_TBI="${CLINVAR_ANNOTATION}.tbi"
GENE_LIST="${DATA_DIR}/pathogenic_genes.txt"
GENE_BED="${DATA_DIR}/pathogenic_genes_1000000bp.bed"
GENE_INTERVAL_LIST="${DATA_DIR}/pathogenic_genes_1000000bp.interval_list"
START_INDEX="${START_INDEX:-0}"
BATCH_SIZE=100

readonly CLINVAR_BASE_PATH CLINVAR_ANNOTATION CLINVAR_ANNOTATION_TBI GENE_LIST GENE_BED GENE_INTERVAL_LIST START_INDEX BATCH_SIZE

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

mapfile -t interval_files < <(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.interval_list 2>/dev/null)
vcf_files=()

if (( ${#interval_files[@]} > 0 )); then
  echo "Found ${#interval_files[@]} interval_list file(s); screening in batches of ${BATCH_SIZE}."
  excluded_intervals=0
  index=0
  while (( index < ${#interval_files[@]} )); do
    batch_paths=("${interval_files[@]:index:BATCH_SIZE}")
    index=$((index + ${#batch_paths[@]}))

    interval_args=()
    tmp_interval_files=()
    shard_ids=()
    for interval_path in "${batch_paths[@]}"; do
      prefix=$(basename "${interval_path}" .interval_list)
      tmp_file=$(mktemp "${TMP_DIR}/${prefix}.XXXXXX.interval_list")
      gsutil -q -u "${GOOGLE_PROJECT}" cp "${interval_path}" "${tmp_file}"
      interval_args+=("-I" "${tmp_file}")
      tmp_interval_files+=("${tmp_file}")
      shard_ids+=("${prefix}")
    done

    output_dir=$(mktemp -d "${TMP_DIR}/iltools_output.XXXXXX")

    gatk IntervalListTools \
      "${interval_args[@]}" \
      -SI "${GENE_INTERVAL_LIST}" \
      --ACTION INTERSECT \
      -O "${output_dir}" \
      --QUIET true

    for shard_id in "${shard_ids[@]}"; do
      out_file="${output_dir}/${shard_id}.interval_list"
      if [[ -f "${out_file}" ]]; then
        if tail -n +2 "${out_file}" | grep -q .; then
          vcf_files+=("${CLINVAR_BASE_PATH}${shard_id}.vcf.bgz")
        else
          excluded_intervals=$((excluded_intervals + 1))
        fi
      else
        excluded_intervals=$((excluded_intervals + 1))
      fi
    done

    rm -rf "${output_dir}"
    for tmp_file in "${tmp_interval_files[@]}"; do
      rm -f "${tmp_file}"
    done
  done
  echo "Selected ${#vcf_files[@]} shard(s) after interval screening (skipped ${excluded_intervals})."
fi

if (( ${#vcf_files[@]} == 0 )); then
  echo "No interval-based selection; defaulting to all VCF shards."
  if ! mapfile -t vcf_files < <(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.vcf.bgz 2>/dev/null); then
    echo "Failed to list VCF shards at ${CLINVAR_BASE_PATH}" >&2
    exit 1
  fi
fi

total=${#vcf_files[@]}
if (( total == 0 )); then
  echo "No VCF shards intersect the gene interval list; nothing to do."
  exit 0
fi

echo "Processing ${total} VCF shard(s)."

if ! [[ ${START_INDEX} =~ ^[0-9]+$ ]]; then
  echo "START_INDEX must be a non-negative integer" >&2
  exit 1
fi

if (( START_INDEX >= total )); then
  echo "START_INDEX (${START_INDEX}) exceeds available shard count (${total})." >&2
  exit 1
fi

last_index=$((total - 1))

processed=${START_INDEX}
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

for (( idx=START_INDEX; idx<total; idx++ )); do
  vcf_file="${vcf_files[idx]}"
  echo "Processing shard ${idx}/${last_index}: ${vcf_file}"
  prefix=$(basename "${vcf_file}" .vcf.bgz)
  local_vcf="${TMP_DIR}/${prefix}.vcf.bgz"
  local_tbi="${local_vcf}.tbi"
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"
  output_tbi="${output_file}.tbi"
  tmp_filtered="${TMP_DIR}/${prefix}.filtered.bcf"
  tmp_annotated="${TMP_DIR}/${prefix}.annot.bcf"

  if [[ -f "${output_file}" && -f "${output_tbi}" ]]; then
    echo "Already processed (files present)."
    processed=$((processed + 1))
    continue
  fi

  copy_start=$(date +%s)
  gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}" "${local_vcf}"
  gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}.tbi" "${local_tbi}"
  copy_end=$(date +%s)
  copy_time_total=$((copy_time_total + (copy_end - copy_start)))
  echo "Copied: ${vcf_file}"

  bcftools view -R "${GENE_BED}" -Ob -o "${tmp_filtered}" "${local_vcf}"

  if ! bcftools view -H "${tmp_filtered}" | grep -q .; then
    echo "No variants after BED filter; writing header-only output for ${prefix}."
    bcftools view -h "${local_vcf}" -Oz -o "${output_file}"
    : > "${output_tbi}"
    rm -f "${local_vcf}" "${local_tbi}" "${tmp_filtered}" "${tmp_filtered}.csi"
    skipped=$((skipped + 1))
    processed=$((processed + 1))
    continue
  fi

  bcftools index -f "${tmp_filtered}"

  bc_start=$(date +%s)
  bcftools annotate \
    -a "${CLINVAR_ANNOTATION}" \
    -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
    -Ob -o "${tmp_annotated}" "${tmp_filtered}"
  bcftools view \
    -i "${FILTER_EXPR}" \
    -Oz -o "${output_file}" "${tmp_annotated}"
  tabix -p vcf "${output_file}"
  bc_end=$(date +%s)
  echo "Filtered: ${vcf_file} to ${output_file}"

  rm -f "${local_vcf}" "${local_tbi}" "${tmp_filtered}" "${tmp_filtered}.csi" "${tmp_annotated}" "${tmp_annotated}.csi"

  filter_time_total=$((filter_time_total + (bc_end - bc_start)))
  processed_with_times=$((processed_with_times + 1))
  processed=$((processed + 1))

done

avg_copy=$(format_avg "$copy_time_total" "$processed_with_times")
avg_filter=$(format_avg "$filter_time_total" "$processed_with_times")

echo "Completed: ${processed}/${total} shard(s); skipped ${skipped}."
echo "Average copy time (processed shards): ${avg_copy}s"
echo "Average filter time (processed shards): ${avg_filter}s"
