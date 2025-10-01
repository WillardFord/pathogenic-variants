#!/usr/bin/env bash
set -euo pipefail

# Annotate all batched VCFs with selected ClinVar INFO tags and index outputs.
# - Input directory:  output/gene_filtered/batched_vcfs
# - Output directory: output/gene_filtered/annotated_vcfs
# - Annotation VCF:   data/clinvar_with_chr.vcf.gz
# - Output naming:    batch_X_annotated.vcf.bgz + .tbi

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

DATA_DIR="${PROJECT_ROOT}/data"
INPUT_DIR="${PROJECT_ROOT}/output/gene_filtered/batched_vcfs"
OUTPUT_DIR="${PROJECT_ROOT}/output/gene_filtered/annotated_vcfs"

ANNOT_VCF="${DATA_DIR}/clinvar_with_chr.vcf.gz"

readonly DATA_DIR INPUT_DIR OUTPUT_DIR ANNOT_VCF

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found in PATH" >&2
  exit 1
fi
if ! command -v tabix >/dev/null 2>&1; then
  echo "tabix not found in PATH" >&2
  exit 1
fi

if [[ ! -f "${ANNOT_VCF}" ]]; then
  echo "Annotation VCF not found: ${ANNOT_VCF}" >&2
  exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Input directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

shopt -s nullglob

# Process both .vcf.bgz and .vcf.gz just in case
mapfile -t inputs < <(find "${INPUT_DIR}" -type f \( -name 'batch_*.vcf.bgz' -o -name 'batch_*.vcf.gz' \) | sort)

if (( ${#inputs[@]} == 0 )); then
  echo "No input VCFs found in ${INPUT_DIR} matching batch_*.vcf.bgz or batch_*.vcf.gz" >&2
  exit 0
fi

echo "Annotating ${#inputs[@]} file(s) from ${INPUT_DIR} -> ${OUTPUT_DIR}"

for in_vcf in "${inputs[@]}"; do
  base=$(basename "${in_vcf}")
  if [[ "${base}" == *.vcf.bgz ]]; then
    out_vcf="${OUTPUT_DIR}/${base%.vcf.bgz}_annotated.vcf.bgz"
  elif [[ "${base}" == *.vcf.gz ]]; then
    out_vcf="${OUTPUT_DIR}/${base%.vcf.gz}_annotated.vcf.bgz"
  else
    # Shouldn't happen given the find filter, but guard anyway
    echo "Skipping unsupported file: ${in_vcf}" >&2
    continue
  fi

  echo "-> ${base}"
  bcftools annotate \
    -a "${ANNOT_VCF}" \
    -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
    -k \
    -Oz -o "${out_vcf}" \
    "${in_vcf}"

  tabix -f -p vcf "${out_vcf}"
done

echo "Done. Outputs in: ${OUTPUT_DIR}"

