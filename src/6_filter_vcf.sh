#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
INPUT_DIR="${PROJECT_ROOT}/output/annotated_vcfs_with_overlapping_ranges"
OUTPUT_DIR="${PROJECT_ROOT}/output/filtered_vcfs_with_overlapping_ranges"
GENE_LIST="${PROJECT_ROOT}/data/pathogenic_genes.txt"
TMP_DIR="${PROJECT_ROOT}/tmp"
JOBS=${JOBS:-2}
BCFTOOLS_THREADS=${BCFTOOLS_THREADS:-2}
CLNSIG_MATCH=${CLNSIG_MATCH:-'INFO/CLNSIG=="Pathogenic,Likely_pathogenic"'}
COMBINED_OUTPUT="${OUTPUT_DIR}/combined_filtered.vcf.bgz"

declare -a TMP_CLEANUP=()
cleanup_all() {
  if (( ${#TMP_CLEANUP[@]} > 0 )); then
    rm -f "${TMP_CLEANUP[@]}" 2>/dev/null || true
  fi
}
trap cleanup_all EXIT

if ! command -v bcftools >/dev/null 2>&1; then
  echo "bcftools not found in PATH." >&2
  exit 1
fi

if ! command -v tabix >/dev/null 2>&1; then
  echo "tabix not found in PATH." >&2
  exit 1
fi

if [[ ! -d "${INPUT_DIR}" ]]; then
  echo "Input directory not found: ${INPUT_DIR}" >&2
  exit 1
fi

if [[ ! -f "${GENE_LIST}" ]]; then
  echo "Gene list not found: ${GENE_LIST}" >&2
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

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"

sanitized_gene_list=$(mktemp "${TMP_DIR}/pathogenic_genes.XXXXXX.txt")
TMP_CLEANUP+=("${sanitized_gene_list}")
# This literally does nothing as the gene list is already sanitized. But it doesn't hurt lols.
grep -Ev '^[[:space:]]*(#|$)' "${GENE_LIST}" > "${sanitized_gene_list}"

if [[ ! -s "${sanitized_gene_list}" ]]; then
  echo "Gene list ${GENE_LIST} is empty after removing blank/comment lines." >&2
  exit 1
fi

readarray -t vcf_files < <(find "${INPUT_DIR}" -maxdepth 1 -type f -name '*.vcf.bgz' | sort)

if (( ${#vcf_files[@]} == 0 )); then
  echo "No annotated VCFs found in ${INPUT_DIR}; nothing to filter." >&2
  exit 0
fi

echo "Filtering ${#vcf_files[@]} file(s) from ${INPUT_DIR} into ${OUTPUT_DIR} using ${JOBS} parallel job(s)..."

n=0
for vcf_path in "${vcf_files[@]}"; do
  base=$(basename "${vcf_path}")
  output_path="${OUTPUT_DIR}/${base}"
  output_tbi="${output_path}.tbi"

  if [[ -f "${output_path}" && -f "${output_tbi}" ]]; then
    echo "Skipping ${base}: filtered output already exists."
    continue
  fi

  {
    header_tmp=$(mktemp "${TMP_DIR}/${base}.header.XXXXXX")
    body_tmp=$(mktemp "${TMP_DIR}/${base}.body.XXXXXX")
    combined_tmp=$(mktemp "${TMP_DIR}/${base}.combined.XXXXXX.vcf")

    cleanup() {
      rm -f "${header_tmp}" "${body_tmp}" "${combined_tmp}" "${combined_tmp}.gz" "${combined_tmp}.gz.tbi"
    }
    trap cleanup EXIT

    bcftools view -h "${vcf_path}" > "${header_tmp}"

    bcftools filter -i "${CLNSIG_MATCH}" --threads "${BCFTOOLS_THREADS}" "${vcf_path}" \
      | bcftools view -H \
      | grep -F -f "${sanitized_gene_list}" > "${body_tmp}" || true

    cat "${header_tmp}" "${body_tmp}" > "${combined_tmp}"

    bcftools view --threads "${BCFTOOLS_THREADS}" -Oz -o "${combined_tmp}.gz" "${combined_tmp}"
    tabix -f -p vcf "${combined_tmp}.gz"

    mv -f "${combined_tmp}.gz" "${output_path}"
    mv -f "${combined_tmp}.gz.tbi" "${output_tbi}"
    cleanup
  } &

  (( ++n % JOBS == 0 )) && wait
done

wait

readarray -t filtered_files < <(find "${OUTPUT_DIR}" -maxdepth 1 -type f -name '*.vcf.bgz' | sort)

if (( ${#filtered_files[@]} == 0 )); then
  echo "No filtered VCFs produced; skipping concatenation." >&2
  exit 0
fi

file_list=$(mktemp "${TMP_DIR}/filtered_file_list.XXXXXX.txt")
TMP_CLEANUP+=("${file_list}")
printf "%s\n" "${filtered_files[@]}" > "${file_list}"

echo "Concatenating ${#filtered_files[@]} filtered VCF(s) into ${COMBINED_OUTPUT}..."
bcftools concat -f "${file_list}" -Oz -o "${COMBINED_OUTPUT}"
tabix -f -p vcf "${COMBINED_OUTPUT}"

echo "Filtering complete. Combined output: ${COMBINED_OUTPUT}"
