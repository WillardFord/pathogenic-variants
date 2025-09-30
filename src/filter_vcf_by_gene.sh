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
GENE_BED="${DATA_DIR}/pathogenic_genes_10000bp.bed"

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

for vcf_file in "${vcf_files[@]}"; do
  prefix=$(basename "${vcf_file}" .vcf.bgz)
  local_vcf="${TMP_DIR}/${prefix}.vcf.bgz"
  local_tbi="${local_vcf}.tbi"
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"
  output_tbi="${output_file}.tbi"

  if [[ -f "${output_file}" && -f "${output_tbi}" ]]; then
    ((processed++))
    ((skipped++))
  else
    copy_start=$(date +%s)
    gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}" "${local_vcf}"
    gsutil -q -u "${GOOGLE_PROJECT}" cp "${vcf_file}.tbi" "${local_tbi}"
    copy_end=$(date +%s)

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

    rm -f "${local_vcf}" "${local_tbi}"

    copy_time_total=$((copy_time_total + (copy_end - copy_start)))
    filter_time_total=$((filter_time_total + (bc_end - bc_start)))
    ((processed_with_times++))
    ((processed++))
  fi

  if (( processed % 10 == 0 )); then
    if (( processed_with_times > 0 )); then
      avg_copy=$(echo "scale=2; ${copy_time_total} / ${processed_with_times}" | bc)
      avg_filter=$(echo "scale=2; ${filter_time_total} / ${processed_with_times}" | bc)
    else
      avg_copy="n/a"
      avg_filter="n/a"
    fi
    echo "Progress: ${processed}/${total} | avg copy: ${avg_copy}s | avg filter: ${avg_filter}s"
  fi

done

if (( processed % 10 != 0 )); then
  if (( processed_with_times > 0 )); then
    avg_copy=$(echo "scale=2; ${copy_time_total} / ${processed_with_times}" | bc)
    avg_filter=$(echo "scale=2; ${filter_time_total} / ${processed_with_times}" | bc)
  else
    avg_copy="n/a"
    avg_filter="n/a"
  fi
  echo "Progress: ${processed}/${total} | avg copy: ${avg_copy}s | avg filter: ${avg_filter}s"
fi

echo "Completed: ${processed}/${total} shard(s); skipped ${skipped}."
