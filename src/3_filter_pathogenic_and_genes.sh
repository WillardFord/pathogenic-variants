#!/usr/bin/env bash
set -euo pipefail

# Filter annotated VCFs to Pathogenic/Likely_pathogenic variants within a gene list.
# - Input directory:  output/gene_filtered/annotaions or output/gene_filtered/annotated_vcfs
# - Output directory: output/gene_filtered/filtered_pathogenic_vcfs
# - Gene list:        data/pathogenic_genes.txt (expects entries with trailing ':', e.g., PMS2:)
# - Criteria:
#     INFO/CLNSIG contains Pathogenic or Likely_pathogenic (token-aware)
#     AND INFO/GENEINFO matches any gene in the list
# - Parallelization: per-file via xargs -P; per-job threads via bcftools --threads

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

DATA_DIR="${PROJECT_ROOT}/data"
DEFAULT_INPUT_A="${PROJECT_ROOT}/output/gene_filtered/annotaions"
DEFAULT_INPUT_B="${PROJECT_ROOT}/output/gene_filtered/annotated_vcfs"
OUTPUT_DIR="${PROJECT_ROOT}/output/gene_filtered/filtered_pathogenic_vcfs"

GENE_LIST="${DATA_DIR}/pathogenic_genes.txt"

# Concurrency controls
DEFAULT_JOBS=$( { command -v nproc >/dev/null 2>&1 && nproc; } || { command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu; } || echo 4 )
JOBS=${JOBS:-${DEFAULT_JOBS}}
BCFTOOLS_THREADS_PER_JOB=${BCFTOOLS_THREADS_PER_JOB:-2}
FORCE=${FORCE:-0}

readonly DATA_DIR DEFAULT_INPUT_A DEFAULT_INPUT_B OUTPUT_DIR GENE_LIST JOBS BCFTOOLS_THREADS_PER_JOB FORCE

command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found in PATH" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found in PATH" >&2; exit 1; }

# Determine input directory (allow override via INPUT_DIR)
INPUT_DIR=${INPUT_DIR:-}
if [[ -z "${INPUT_DIR}" ]]; then
  if [[ -d "${DEFAULT_INPUT_A}" ]]; then
    INPUT_DIR="${DEFAULT_INPUT_A}"
  elif [[ -d "${DEFAULT_INPUT_B}" ]]; then
    INPUT_DIR="${DEFAULT_INPUT_B}"
  else
    echo "No input directory found. Checked: ${DEFAULT_INPUT_A} and ${DEFAULT_INPUT_B}. Set INPUT_DIR explicitly." >&2
    exit 1
  fi
fi

[[ -f "${GENE_LIST}" ]] || { echo "Gene list not found: ${GENE_LIST}" >&2; exit 1; }

mkdir -p "${OUTPUT_DIR}"

# Build gene regex: join lines into G1:|G2:|... (list already contains trailing ':')
GENE_PATTERN=$(awk 'NF {print}' "${GENE_LIST}" | paste -sd'|' -)
if [[ -z "${GENE_PATTERN}" ]]; then
  echo "Gene list produced an empty pattern" >&2
  exit 1
fi

# Token-aware CLNSIG match: ensure Pathogenic or Likely_pathogenic matches as a whole token within '|' separated values
CLNSIG_EXPR='INFO/CLNSIG ~ "(^|\\|)(Pathogenic|Likely_pathogenic)($|\\|)"'
GENEINFO_EXPR="INFO/GENEINFO ~ \"(${GENE_PATTERN})\""
BCF_EXPR="${CLNSIG_EXPR} && ${GENEINFO_EXPR}"

# Count inputs
NUM_INPUTS=$(find "${INPUT_DIR}" -type f \( -name 'batch_*_annotated.vcf.bgz' -o -name 'batch_*.vcf.bgz' -o -name '*.vcf.bgz' \) | wc -l | awk '{print $1}')
if (( NUM_INPUTS == 0 )); then
  echo "No input VCFs found in ${INPUT_DIR}" >&2
  exit 0
fi

echo "Filtering ${NUM_INPUTS} file(s) from ${INPUT_DIR} -> ${OUTPUT_DIR} using ${JOBS} job(s); bcftools --threads=${BCFTOOLS_THREADS_PER_JOB}."
echo "bcftools -i expression: ${BCF_EXPR}"

export OUTPUT_DIR BCF_EXPR BCFTOOLS_THREADS_PER_JOB FORCE

find "${INPUT_DIR}" -type f \( -name 'batch_*_annotated.vcf.bgz' -o -name 'batch_*.vcf.bgz' -o -name '*.vcf.bgz' \) -print0 \
  | xargs -0 -n 1 -P "${JOBS}" bash -c '
      set -euo pipefail
      in_vcf="$1"
      base=$(basename "$in_vcf")
      out_vcf="${OUTPUT_DIR}/${base%*.vcf.bgz}_pathogenic.vcf.bgz"

      if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" && "${FORCE}" != "1" ]]; then
        echo "[SKIP] ${base} -> already filtered"
        exit 0
      fi

      echo "[FILTER] ${base} -> $(basename "$out_vcf")"
      tmp_out="${out_vcf}.tmp"
      trap "rm -f \"$tmp_out\" \"${tmp_out}.tbi\"" EXIT

      bcftools view \
        -i "${BCF_EXPR}" \
        --threads "${BCFTOOLS_THREADS_PER_JOB}" \
        -Oz -o "$tmp_out" \
        "$in_vcf"

      tabix -f -p vcf "$tmp_out"
      mv -f "$tmp_out" "$out_vcf"
      mv -f "${tmp_out}.tbi" "${out_vcf}.tbi"
    ' _

echo "Done. Outputs in: ${OUTPUT_DIR}"

