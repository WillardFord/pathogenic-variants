# For some reason the scripts don't work with BCF tool filtering.
# So just copy and paste below:
:'
JOBS=${JOBS:-16}
BCFTOOLS_THREADS_PER_JOB=${BCFTOOLS_THREADS_PER_JOB:-2}
INPUT_DIR="output/annotated_vcfs_with_overlapping_ranges"
OUTPUT_DIR="output/filtered_vcfs_with_overlapping_ranges"

mkdir -p "${OUTPUT_DIR}"

n=0
for f in "${INPUT_DIR}"/*.vcf.bgz; do
  base="$(basename "$f")"
  out="${OUTPUT_DIR}/${base}"
  {
    bcftools view -h "$f" > "${out}.tmp"
    bcftools filter -i 'INFO/CLNSIG=="Pathogenic,Likely_pathogenic"'  "$f" | bcftools view -H | grep -Ff data/pathogenic_genes.txt >> "${out}.tmp"
    bcftools view -Oz -o "${out}" "${out}.tmp"
    tabix -f -p vcf "${out}"
    rm -f "${out}.tmp"
  } &
  (( ++n % JOBS == 0 )) && wait
done
wait
'

"
# Concatonation
OUTPUT_DIR="output/"
ls -1v "${OUTPUT_DIR}"/*.vcf.bgz > "${OUTPUT_DIR}/file_list.txt"
bcftools concat -f "${OUTPUT_DIR}/file_list.txt" -Oz -o "${OUTPUT_DIR}/combined.vcf.gz"
bcftools index "${OUTPUT_DIR}/combined.vcf.gz"
"