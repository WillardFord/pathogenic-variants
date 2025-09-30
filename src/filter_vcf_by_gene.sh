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
START_INDEX="${START_INDEX:-0}"

readonly CLINVAR_BASE_PATH CLINVAR_ANNOTATION CLINVAR_ANNOTATION_TBI GENE_LIST GENE_BED START_INDEX

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

echo "START_INDEX set to ${START_INDEX}"

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

if (( START_INDEX >= total )); then
  echo "START_INDEX (${START_INDEX}) exceeds available shards (${total})." >&2
  exit 1
fi

processed=$START_INDEX
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

for (( idx=START_INDEX; idx<total; idx++ )); do
  vcf_file="${vcf_files[idx]}"
  echo "Processing shard ${idx}/${total - 1}: ${vcf_file}"
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
  echo "Copied: ${vcf_file}"

  bcftools view -R "${GENE_BED}" -Ob -o "${tmp_filtered}" "${local_vcf}"

  if ! bcftools view -H "${tmp_filtered}" | grep -q .; then
    echo "No variants after BED filter; creating empty marker."
    rm -f "${local_vcf}" "${local_tbi}" "${tmp_filtered}" "${tmp_filtered}.csi"
    : > "${output_file}"  # create empty zero-byte marker
    : > "${output_tbi}"
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

  copy_time_total=$((copy_time_total + (copy_end - copy_start)))
  filter_time_total=$((filter_time_total + (bc_end - bc_start)))
  processed_with_times=$((processed_with_times + 1))
  processed=$((processed + 1))

done

avg_copy=$(format_avg "$copy_time_total" "$processed_with_times")
avg_filter=$(format_avg "$filter_time_total" "$processed_with_times")

echo "Completed: ${processed}/${total} shard(s); skipped ${skipped}."
echo "Average copy time: ${avg_copy}s"
echo "Average filter time: ${avg_filter}s"
