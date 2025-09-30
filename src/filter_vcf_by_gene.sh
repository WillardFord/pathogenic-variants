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

readonly CLINVAR_BASE_PATH CLINVAR_ANNOTATION CLINVAR_ANNOTATION_TBI GENE_LIST GENE_BED

for path in "${CLINVAR_ANNOTATION}" "${CLINVAR_ANNOTATION_TBI}" "${GENE_LIST}" "${GENE_BED}"; do
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
echo "Using region BED: ${GENE_BED}"

vcf_files=()
while IFS= read -r line; do
  [[ -n "${line}" ]] && vcf_files+=("${line}")
done < <(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.vcf.bgz)

total=${#vcf_files[@]}
if (( total == 0 )); then
  echo "No VCF shards found at ${CLINVAR_BASE_PATH}" >&2
  exit 0
fi

echo "Discovered ${total} VCF shard(s) to evaluate."

processed=0
processed_with_times=0
skipped=0
copy_time_total=0
filter_time_total=0

format_avg() {
  local total=$1
  local count=$2
  if (( count == 0 )); then
    printf 'n/a'
  else
    awk -v t="$total" -v c="$count" 'BEGIN { printf "%.2f", t / c }'
  fi
}

report_progress() {
  local label=$1
  local avg_copy=$(format_avg "$copy_time_total" "$processed_with_times")
  local avg_filter=$(format_avg "$filter_time_total" "$processed_with_times")
  echo "Progress: ${processed}/${total} | ${label} | avg copy: ${avg_copy}s | avg filter: ${avg_filter}s"
}

for vcf_file in "${vcf_files[@]}"; do
  echo "Processing: ${vcf_file}"
  prefix=$(basename "${vcf_file}" .vcf.bgz)
  local_vcf="${TMP_DIR}/${prefix}.vcf.bgz"
  local_tbi="${local_vcf}.tbi"
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"
  output_tbi="${output_file}.tbi"

  if [[ -f "${output_file}" && -f "${output_tbi}" ]]; then
    ((processed++))
    ((skipped++))
    if (( processed % 10 == 0 )); then
      report_progress "skipping"
    fi
    continue
  fi

  copy_start=$(date +%s)
  gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}" "${local_vcf}"
  gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}.tbi" "${local_tbi}"
  copy_end=$(date +%s)
  echo "Copied: ${vcf_file} to ${local_vcf}"
  bc_start=$(date +%s)
  bcftools view -R "${GENE_BED}" -Ou "${local_vcf}" \
    | bcftools annotate \
        -a "${CLINVAR_ANNOTATION}" \
        -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
        -Ou - \
    | bcftools view \
        -i "${FILTER_EXPR}" \
        -Oz -o "${output_file}"
  tabix -p vcf "${output_file}"
  bc_end=$(date +%s)
  echo "Filtered: ${vcf_file} to ${output_file}"
  rm -f "${local_vcf}" "${local_tbi}"

  copy_time_total=$((copy_time_total + (copy_end - copy_start)))
  filter_time_total=$((filter_time_total + (bc_end - bc_start)))
  ((processed_with_times++))
  ((processed++))

  if (( processed % 10 == 0 )); then
    report_progress "processing"
  fi

done

if (( processed % 10 != 0 )); then
  report_progress "final"
fi

echo "Completed: ${processed}/${total} shard(s); skipped ${skipped}."
