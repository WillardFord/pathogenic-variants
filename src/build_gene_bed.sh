#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DATA_DIR="${PROJECT_ROOT}/data"

GENE_LIST_DEFAULT="${DATA_DIR}/pathogenic_genes.txt"
GTF_DEFAULT="${DATA_DIR}/gencode.v49.primary_assembly.annotation.gtf.gz"
FLANK_DEFAULT=1000000
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

if ! command -v bedtools >/dev/null 2>&1; then
  echo "bedtools not found in PATH" >&2
  exit 1
fi

TMP_DIR=$(mktemp -d)
trap 'rm -rf "${TMP_DIR}"' EXIT

GENE_NAMES="${TMP_DIR}/genes.txt"
RAW_BED="${TMP_DIR}/raw.bed"

sed 's/:$//' "${GENE_LIST}" > "${GENE_NAMES}"

zcat -f "${GTF_PATH}" \
  | awk -F'\t' -v OFS='\t' -v flank="${FLANK}" -v list="${GENE_NAMES}" '
      BEGIN {
        while ((getline name < list) > 0) {
          gsub(/^[ \t]+|[ \t]+$/, "", name);
          if (name != "") keep[name] = 1;
        }
      }
      $3 == "gene" {
        gene = "";
        split($9, attrs, ";");
        for (i = 1; i <= length(attrs); i++) {
          entry = attrs[i];
          gsub(/^[ \t]+|[ \t]+$/, "", entry);
          if (index(entry, "gene_name") == 1) {
            split(entry, parts, "\"");
            if (length(parts) >= 2) {
              gene = parts[2];
              break;
            }
          }
        }
        if (gene != "" && (gene in keep)) {
          start = $4 - 1 - flank;
          if (start < 0) start = 0;
          end = $5 + flank;
          print $1, start, end, gene;
        }
      }
    ' > "${RAW_BED}"

if [[ ! -s "${RAW_BED}" ]]; then
  echo "No gene intervals found for the provided list" >&2
  exit 1
fi

bedtools sort -i "${RAW_BED}" \
  | bedtools merge -i - -c 4 -o distinct \
  > "${OUTPUT_PATH}"

echo "Wrote BED to ${OUTPUT_PATH}"
