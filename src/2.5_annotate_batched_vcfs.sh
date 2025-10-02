#!/usr/bin/env bash
#set -euo pipefail

# Parallel annotation of all batched VCFs with selected ClinVar INFO tags.
# - Input:  output/gene_filtered/batched_vcfs/batch_*.vcf.bgz(.gz)
# - Output: output/gene_filtered/annotated_vcfs/batch_X_annotated.vcf.bgz + .tbi
# Controls:
#   JOBS                     .. number of parallel files (default: CPU cores)
#   BCFTOOLS_THREADS_PER_JOB .. threads passed to bcftools per file (default: 2)
#   FORCE=1                  .. re-run even if outputs exist


# AHHHHH this doesn't work. But when I copy paste the below command into the terminal
# It works fine. What is the difference?
# bcftools annotate -a data/clinvar_with_chr.vcf.gz -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC -k -Ou output/gene_filtered/batched_vcfs/batch_10.vcf.bgz
:'Okay I gave up and just used the below command in terminal and it worked.
JOBS=8
n=0
for f in output/gene_filtered/batched_vcfs/*.vcf.bgz; do
  base=$(basename "$f")
  out="output/gene_filtered/annotated_vcfs/$base"
  {
  bcftools annotate -a data/clinvar_with_chr.vcf.gz \
    -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
    -k --threads 2 -Oz -o "$out" "$f" \
    && tabix -f -p vcf "$out"
  } &
  (( ++n % JOBS == 0 )) && wait
done
wait
# This will create defunct process if you are not careful about how you kill them
# So just let them finish to be sure. Shouldnt take too long.
'

:'
SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

DATA_DIR="${PROJECT_ROOT}/data"
INPUT_DIR="${PROJECT_ROOT}/output/gene_filtered/batched_vcfs"
OUTPUT_DIR="${PROJECT_ROOT}/output/gene_filtered/annotated_vcfs"

ANNOT_VCF="${DATA_DIR}/clinvar_with_chr.vcf.gz"

# Derive sensible defaults for parallelism
DEFAULT_JOBS=$( { command -v nproc >/dev/null 2>&1 && nproc; } || { command -v sysctl >/dev/null 2>&1 && sysctl -n hw.ncpu; } || echo 4 )
JOBS=${JOBS:-${DEFAULT_JOBS}}
BCFTOOLS_THREADS_PER_JOB=${BCFTOOLS_THREADS_PER_JOB:-2}
FORCE=${FORCE:-0}

readonly DATA_DIR INPUT_DIR OUTPUT_DIR ANNOT_VCF JOBS BCFTOOLS_THREADS_PER_JOB FORCE

# Tool checks
command -v bcftools >/dev/null 2>&1 || { echo "bcftools not found in PATH" >&2; exit 1; }
command -v tabix >/dev/null 2>&1 || { echo "tabix not found in PATH" >&2; exit 1; }

# Inputs sanity
[[ -f "${ANNOT_VCF}" ]] || { echo "Annotation VCF not found: ${ANNOT_VCF}" >&2; exit 1; }
if [[ ! -f "${ANNOT_VCF}.tbi" && ! -f "${ANNOT_VCF}.csi" ]]; then
  echo "Annotation VCF index not found: ${ANNOT_VCF}.tbi (or .csi)" >&2
  exit 1
fi
[[ -d "${INPUT_DIR}" ]] || { echo "Input directory not found: ${INPUT_DIR}" >&2; exit 1; }
mkdir -p "${OUTPUT_DIR}"

# Count inputs
#NUM_INPUTS=$(find "${INPUT_DIR}" -type f \( -name 'batch_*.vcf.bgz' -o -name 'batch_*.vcf.gz' \) | wc -l | awk '{print $1}')
'
:'if (( NUM_INPUTS == 0 )); then
  echo "No input VCFs found in ${INPUT_DIR} matching batch_*.vcf.bgz or batch_*.vcf.gz" >&2
  exit 0
fi

echo "Annotating ${NUM_INPUTS} file(s) from ${INPUT_DIR} -> ${OUTPUT_DIR} using ${JOBS} parallel job(s); bcftools --threads=${BCFTOOLS_THREADS_PER_JOB}."

export OUTPUT_DIR ANNOT_VCF BCFTOOLS_THREADS_PER_JOB FORCE

# Parallel per-file processing via xargs -P
'
#find "${INPUT_DIR}" -type f \( -name 'batch_*.vcf.bgz' -o -name 'batch_*.vcf.gz' \) -print0 \
#  | xargs -0 -n 1 -P "${JOBS}" bash -c '
#      set -euo pipefail
#      in_vcf="$1"
#      base=$(basename "$in_vcf")
#      if [[ "$base" == *.vcf.bgz ]]; then
#        out_vcf="${OUTPUT_DIR}/${base%.vcf.bgz}_annotated.vcf.bgz"
#      elif [[ "$base" == *.vcf.gz ]]; then
#        out_vcf="${OUTPUT_DIR}/${base%.vcf.gz}_annotated.vcf.bgz"
#      else
#        echo "[SKIP] Unsupported file: $in_vcf" >&2
#        exit 0
#      fi

#      if [[ -s "$out_vcf" && -s "${out_vcf}.tbi" && "${FORCE}" != "1" ]]; then
#        echo "[SKIP] ${base} -> already annotated"
#        exit 0
#      fi
#
#      echo "[ANNOTATE] ${base} -> $(basename "$out_vcf")"
#      tmp_out="${out_vcf}.tmp"
#      trap "rm -f \"$tmp_out\" \"${tmp_out}.tbi\"" EXIT
#
#      bcftools annotate \
#        -a "${ANNOT_VCF}" \
#        -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
#        -k \
#        --threads "${BCFTOOLS_THREADS_PER_JOB}" \
#        -Oz -o "$tmp_out" \
#        "$in_vcf"
#
#      tabix -f -p vcf "$tmp_out"
#      mv -f "$tmp_out" "$out_vcf"
#      mv -f "${tmp_out}.tbi" "${out_vcf}.tbi"
#    ' _
#
#echo "Done. Outputs in: ${OUTPUT_DIR}"
#'


