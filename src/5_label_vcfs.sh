#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
INPUT_DIR="${PROJECT_ROOT}/output/downloaded_vcfs_with_overlapping_ranges"
OUTPUT_DIR="${PROJECT_ROOT}/output/annotated_vcfs_with_overlapping_ranges"
ANNOTATION="${PROJECT_ROOT}/data/clinvar_with_chr.vcf.gz"
ANNOTATION_TBI="${ANNOTATION}.tbi"
JOBS=${JOBS:-16}
BCFTOOLS_THREADS=${BCFTOOLS_THREADS:-2}

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found in PATH." >&2
  exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
  echo "tabix not found in PATH." >&2
  exit 1
fi

if [[ ! -f "${ANNOTATION}" ]]; then
  echo "ClinVar annotation file not found: ${ANNOTATION}" >&2
  exit 1
fi

if [[ ! -f "${ANNOTATION_TBI}" ]]; then
  echo "ClinVar annotation index not found: ${ANNOTATION_TBI}" >&2
  exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Input directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

if ! [[ "${JOBS}" =~ ^[0-9]+$ ]] || (( JOBS <= 0 )); then
  echo "JOBS must be a positive integer (got: ${JOBS})." >&2
  exit 1
fi

if ! [[ "${BCFTOOLS_THREADS}" =~ ^[0-9]+$ ]] || (( BCFTOOLS_THREADS <= 0 )); then
  echo "BCFTOOLS_THREADS must be a positive integer (got: ${BCFTOOLS_THREADS})." >&2
  exit 1
fi

mkdir -p "${OUTPUT_DIR}"

readarray -t vcf_files < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.vcf.bgz' -print | sort)

if (( ${#vcf_files[@]} == 0 )); then
  echo "No VCF files (*.vcf.bgz) found in ${INPUT_DIR}." >&2
  exit 0
fi

echo "Annotating ${#vcf_files[@]} file(s) from ${INPUT_DIR} into ${OUTPUT_DIR} using ${JOBS} parallel job(s)..."

n=0
for vcf_path in "${vcf_files[@]}"; do
  base=$(basename "${vcf_path}")
  output_path="${OUTPUT_DIR}/${base%.vcf.bgz}_annotated.vcf.bgz"
  output_tbi="${output_path}.tbi"

  if [[ -f "${output_path}" && -f "${output_tbi}" ]]; then
    echo "Skipping ${base}: annotated output already exists."
    continue
  fi

  {
    bcftools annotate \
      -a "${ANNOTATION}" \
      -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
      -k \
      --threads "${BCFTOOLS_THREADS}" \
      -Oz \
      -o "${output_path}" \
      "${vcf_path}"

    tabix -f -p vcf "${output_path}"
  } &

  (( ++n % JOBS == 0 )) && wait
done

wait

echo "Annotation complete."
