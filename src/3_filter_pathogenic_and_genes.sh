# For some reason the scripts don't work with BCF tool filtering.
# So just copy and paste below:
:'
JOBS=${JOBS:-8}
BCFTOOLS_THREADS_PER_JOB=${BCFTOOLS_THREADS_PER_JOB:-2}
INPUT_DIR="output/gene_filtered/annotated_vcfs/"
OUTPUT_DIR="output/gene_filtered/filtered_pathogenic_vcfs"

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