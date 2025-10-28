#!/bin/bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PROJECT_ROOT=$(cd "${SCRIPT_DIR}/.." && pwd)
DOWNLOAD_LIST="${PROJECT_ROOT}/output/vcf_download_list.txt"
DEST_DIR="${PROJECT_ROOT}/output/downloaded_vcfs_with_overlapping_ranges"

mkdir -p "${DEST_DIR}"

gsutil -u "${GOOGLE_PROJECT}" -m \
  -o "GSUtil:parallel_thread_count=16" \
  -o "GSUtil:parallel_process_count=1" \
  cp -I "${DEST_DIR}/" < "${DOWNLOAD_LIST}"
