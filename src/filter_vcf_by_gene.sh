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

for path in "${CLINVAR_ANNOTATION}" "${CLINVAR_ANNOTATION_TBI}" "${GENE_LIST}"; do
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

vcf_files=$(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.vcf.bgz)

for vcf_file in ${vcf_files}; do
  prefix=$(basename "${vcf_file}" .vcf.bgz)
  local_vcf="${TMP_DIR}/${prefix}.vcf.bgz"
  local_tbi="${local_vcf}.tbi"
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"

  echo "Processing: ${vcf_file}"
  echo "Output: ${output_file}"

  gsutil -u "${GOOGLE_PROJECT}" cp "${vcf_file}" "${local_vcf}"
  gsutil -u "${GOOGLE_PROJECT}" cp "${vcf_file}.tbi" "${local_tbi}"

  bcftools annotate \
    -a "${CLINVAR_ANNOTATION}" \
    -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
    -Ou "${local_vcf}" \
  | bcftools view \
    -i "${FILTER_EXPR}" \
    -Oz -o "${output_file}"

  tabix -p vcf "${output_file}"

  rm -f "${local_vcf}" "${local_tbi}"

  echo "Completed: ${output_file}"
  echo "---"

  break

done

echo "All VCF files processed successfully!"
