#!/usr/bin/env bash
set -euo pipefail

# This script screens every interval_list shard in output/interval_lists
# for overlap against the configured gene interval list. Shards that overlap
# have their matching remote VCF paths recorded in output/vcf_download_list.txt,
# and their individual overlap interval_lists are written to
# data/interval_lists_with_overlapping_ranges/.

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)

DATA_DIR="${PROJECT_ROOT}/data"
OUTPUT_DIR="${PROJECT_ROOT}/output"
TMP_DIR="${PROJECT_ROOT}/tmp"

INTERVAL_LIST_DIR="${OUTPUT_DIR}/interval_lists"
INTERMEDIATE_DIR="${DATA_DIR}/interval_lists_with_overlapping_ranges"
DOWNLOAD_LIST="${OUTPUT_DIR}/vcf_download_list.txt"
INTERVAL_LISTS_WITH_OVERLAP="${OUTPUT_DIR}/interval_lists_with_overlaps.txt"
GENE_INTERVAL_LIST="${DATA_DIR}/pathogenic_genes_1000000bp.interval_list"
REMOTE_VCF_BASE="${REMOTE_VCF_BASE:-gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf}"
MAX_JOBS="${MAX_JOBS:-}"

RESULTS_DIR=""
fifo_fd_open=0

cleanup() {
  if (( fifo_fd_open )); then
    exec 9>&-
    exec 9<&-
    fifo_fd_open=0
  fi
  if [[ -n "${RESULTS_DIR}" && -d "${RESULTS_DIR}" ]]; then
    rm -rf "${RESULTS_DIR}"
  fi
}
trap cleanup EXIT

die() {
  echo "Error: $*" >&2
  exit 1
}

detect_default_jobs() {
  if command -v nproc >/dev/null 2>&1; then
    nproc
    return
  fi
  local uname_out
  uname_out=$(uname -s 2>/dev/null || echo "")
  case "${uname_out}" in
    Darwin)
      sysctl -n hw.ncpu 2>/dev/null || echo 4
      ;;
    *)
      getconf _NPROCESSORS_ONLN 2>/dev/null || echo 4
      ;;
  esac
}

run_overlap() {
  local interval_path="$1"
  local result_dir="$2"
  local base prefix overlap_path hit_file

  base=$(basename "${interval_path}")
  prefix="${base%.interval_list}"
  overlap_path="${INTERMEDIATE_DIR}/${prefix}.overlap.interval_list"
  hit_file="${result_dir}/${prefix}.hit"

  rm -f "${overlap_path}" "${hit_file}"

  gatk IntervalListTools \
    -I "${interval_path}" \
    -SI "${GENE_INTERVAL_LIST}" \
    --ACTION OVERLAPS \
    -O "${overlap_path}" \
    --QUIET true

  if LC_ALL=C grep -q -m1 '^[^@]' "${overlap_path}"; then
    printf '%s\n' "${interval_path}" > "${hit_file}"
  fi
}

[[ -f "${GENE_INTERVAL_LIST}" ]] || die "Missing gene interval list at ${GENE_INTERVAL_LIST}"
[[ -d "${INTERVAL_LIST_DIR}" ]] || die "Missing interval list directory at ${INTERVAL_LIST_DIR}"
command -v gatk >/dev/null 2>&1 || die "gatk not found in PATH"

if [[ -z "${MAX_JOBS}" ]]; then
  MAX_JOBS=$(detect_default_jobs)
fi
case "${MAX_JOBS}" in
  ''|*[!0-9]*)
    die "MAX_JOBS must be a positive integer (received '${MAX_JOBS}')"
    ;;
esac
if (( MAX_JOBS < 1 )); then
  MAX_JOBS=1
fi

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"

if [[ -d "${INTERMEDIATE_DIR}" ]]; then
  find "${INTERMEDIATE_DIR}" -type f -name '*.interval_list' -delete
else
  mkdir -p "${INTERMEDIATE_DIR}"
fi

> "${DOWNLOAD_LIST}"
> "${INTERVAL_LISTS_WITH_OVERLAP}"

interval_files=()
while IFS= read -r path; do
  interval_files+=("${path}")
done < <(find "${INTERVAL_LIST_DIR}" -type f -name '*.interval_list' -print | sort)

if (( ${#interval_files[@]} == 0 )); then
  echo "No interval_list files found under ${INTERVAL_LIST_DIR}; nothing to do."
  exit 0
fi

if (( MAX_JOBS > ${#interval_files[@]} )); then
  MAX_JOBS=${#interval_files[@]}
fi

RESULTS_DIR=$(mktemp -d "${TMP_DIR}/interval_overlap.XXXXXX")

fifo_path=$(mktemp "${TMP_DIR}/interval_sem.XXXXXX")
rm -f "${fifo_path}"
mkfifo "${fifo_path}"
exec 9<>"${fifo_path}"
rm -f "${fifo_path}"
fifo_fd_open=1

job=0
total=${#interval_files[@]}
pids=()

while (( job < MAX_JOBS )); do
  printf 'token\n' >&9
  job=$((job + 1))
done

index=0
for interval_path in "${interval_files[@]}"; do
  read -r _ <&9
  index=$((index + 1))
  echo "[${index}/${total}] Scheduling $(basename "${interval_path}")"
  {
    run_overlap "${interval_path}" "${RESULTS_DIR}"
    printf 'token\n' >&9
  } &
  pids+=("$!")
done

set +e
failures=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    failures=$((failures + 1))
  fi
done
set -e

exec 9>&-
exec 9<&-
fifo_fd_open=0

if (( failures > 0 )); then
  die "${failures} overlap job(s) failed"
fi

overlapping_paths=()
while IFS= read -r -d '' hit_file; do
  while IFS= read -r path; do
    if [[ -n "${path}" ]]; then
      overlapping_paths+=("${path}")
    fi
  done < "${hit_file}"
done < <(find "${RESULTS_DIR}" -type f -name '*.hit' -print0)

if (( ${#overlapping_paths[@]} == 0 )); then
  echo "No overlaps detected across ${total} interval list(s)."
  exit 0
fi

unique_tmp=$(mktemp "${TMP_DIR}/interval_hits_unique.XXXXXX")
printf '%s\n' "${overlapping_paths[@]}" | sort -u > "${unique_tmp}"

unique_paths=()
while IFS= read -r path; do
  if [[ -n "${path}" ]]; then
    unique_paths+=("${path}")
  fi
done < "${unique_tmp}"
rm -f "${unique_tmp}"

remote_base="${REMOTE_VCF_BASE%/}"
[[ -n "${remote_base}" ]] || die "REMOTE_VCF_BASE must not be empty"

> "${DOWNLOAD_LIST}"
> "${INTERVAL_LISTS_WITH_OVERLAP}"

for path in "${unique_paths[@]}"; do
  base=$(basename "${path}")
  prefix="${base%.interval_list}"
  printf '%s/%s.vcf.bgz\n' "${remote_base}" "${prefix}" >> "${DOWNLOAD_LIST}"
  printf '%s/%s.vcf.bgz.tbi\n' "${remote_base}" "${prefix}" >> "${DOWNLOAD_LIST}"
  printf '%s\n' "${path}" >> "${INTERVAL_LISTS_WITH_OVERLAP}"
done

echo "Identified ${#unique_paths[@]} interval list(s) with overlaps."
echo "Wrote VCF download list to ${DOWNLOAD_LIST}."
echo "Wrote matching interval list paths to ${INTERVAL_LISTS_WITH_OVERLAP}."
