#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DATA_DIR="${PROJECT_ROOT}/data"
OUTPUT_DIR="${PROJECT_ROOT}/output"

CLINVAR_BASE_PATH="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf/"
CLINVAR_ANNOTATION="${DATA_DIR}/clinvar.vcf.gz"
GENE_LIST="${DATA_DIR}/pathogenic_genes.txt"

if [[ ! -f "${CLINVAR_ANNOTATION}" ]]; then
  echo "Missing ClinVar annotation VCF at ${CLINVAR_ANNOTATION}" >&2
  exit 1
fi

if [[ ! -f "${GENE_LIST}" ]]; then
  echo "Missing gene list at ${GENE_LIST}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

vcf_files=$(gsutil -u "${GOOGLE_PROJECT}" ls "${CLINVAR_BASE_PATH}"*.vcf.bgz)

for vcf_file in ${vcf_files}; do
  prefix=$(basename "${vcf_file}" .vcf.bgz)
  output_file="${OUTPUT_DIR}/${prefix}_filtered.vcf.bgz"

  echo "Processing: ${vcf_file}"
  echo "Output: ${output_file}"

  gsutil -u "${GOOGLE_PROJECT}" cat "${vcf_file}" \
    | bcftools annotate \
        -a "${CLINVAR_ANNOTATION}" \
        -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
        -Ou - \
    | bcftools view \
        -i "INFO/GENEINFO ~ @${GENE_LIST}" \
        -Oz -o "${output_file}"

  tabix -p vcf "${output_file}"

  echo "Completed: ${output_file}"
  echo "---"

  break

done

echo "All VCF files processed successfully!"
