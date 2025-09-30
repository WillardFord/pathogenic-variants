#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DATA_DIR="${PROJECT_ROOT}/data"

GENE_LIST_DEFAULT="${DATA_DIR}/pathogenic_genes.txt"
GTF_DEFAULT="${DATA_DIR}/gencode.v39.annotation.gtf.gz"
FLANK_DEFAULT=10000
OUTPUT_DEFAULT="${DATA_DIR}/pathogenic_genes_${FLANK_DEFAULT}bp.bed"

usage() {
  cat <<USAGE
Usage: $(basename "$0") [-g gtf.gz] [-l gene_list] [-o output.bed] [-f flank_bp]
  -g  Path to GENCODE (or compatible) GTF file (default: ${GTF_DEFAULT})
  -l  Gene list with symbols per line (default: ${GENE_LIST_DEFAULT})
  -o  Output BED path (default: ${OUTPUT_DEFAULT})
  -f  Flank size in bp to extend on each side (default: ${FLANK_DEFAULT})
USAGE
}

GTF_PATH="${GTF_DEFAULT}"
GENE_LIST="${GENE_LIST_DEFAULT}"
OUTPUT_PATH="${OUTPUT_DEFAULT}"
FLANK="${FLANK_DEFAULT}"

while getopts ":g:l:o:f:h" opt; do
  case ${opt} in
    g) GTF_PATH="${OPTARG}" ;;
    l) GENE_LIST="${OPTARG}" ;;
    o) OUTPUT_PATH="${OPTARG}" ;;
    f) FLANK="${OPTARG}" ;;
    h) usage; exit 0 ;;
    :) echo "Option -${OPTARG} requires an argument" >&2; usage; exit 1 ;;
    \?) echo "Invalid option -${OPTARG}" >&2; usage; exit 1 ;;
  esac
done

if [[ ! -f "${GENE_LIST}" ]]; then
  echo "Missing gene list: ${GENE_LIST}" >&2
  exit 1
fi

if [[ ! -f "${GTF_PATH}" ]]; then
  echo "Missing GTF file: ${GTF_PATH}" >&2
  exit 1
fi

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

GENE_NAMES="${TMP_DIR}/genes.txt"
RAW_BED="${TMP_DIR}/raw.bed"

# Normalise gene names by stripping trailing colon if present
sed 's/:$//' "${GENE_LIST}" > "${GENE_NAMES}"

# Extract gene features from the GTF and pad using awk
zcat -f "${GTF_PATH}" \
  | awk -F'\t' -v OFS='\t' -v flank="${FLANK}" -v genes="${GENE_NAMES}" '
      BEGIN { while ((getline g < genes) > 0) keep[g] = 1 }
      $3 == "gene" {
        if (match($9, /gene_name "([^"]+)"/, m)) {
          gene = m[1]
          if (keep[gene]) {
            start = $4 - 1 - flank
            if (start < 0) start = 0
            end = $5 + flank
            print $1, start, end, gene
          }
        }
      }
    ' > "${RAW_BED}"

if [[ ! -s "${RAW_BED}" ]]; then
  echo "No gene intervals found for the provided list" >&2
  exit 1
fi

# Sort, merge overlapping segments, and keep distinct gene labels
bedtools sort -i "${RAW_BED}" \
  | bedtools merge -i - -c 4 -o distinct \
  > "${OUTPUT_PATH}"

echo "Wrote BED to ${OUTPUT_PATH}"
