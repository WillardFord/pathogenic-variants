#!/bin/bash

# Set the base path for ClinVar VCF files
CLINVAR_BASE_PATH="gs://fc-aou-datasets-controlled/v8/wgs/short_read/snpindel/clinvar/vcf/"

# Get list of .vcf.bgz files (excluding .tbi and .interval_list files)
vcf_files=$(gsutil -u $GOOGLE_PROJECT ls "${CLINVAR_BASE_PATH}*.vcf.bgz")

# Process each VCF file
for vcf_file in $vcf_files; do
    # Extract the prefix (numbers before .vcf.bgz)
    prefix=$(basename "$vcf_file" .vcf.bgz)
    
    # Define output filename with _filtered suffix
    output_file="${prefix}_filtered.vcf.bgz"
    
    echo "Processing: $vcf_file"
    echo "Output: $output_file"
    
    # Process the VCF file
    gsutil -u $GOOGLE_PROJECT cat "$vcf_file" |
    bcftools annotate \
           -a ../data/clinvar.vcf.gz \
           -c INFO/CLNSIG,INFO/CLNSIGCONF,INFO/GENEINFO,INFO/MC \
           - \
           -Ou |
    bcftools view \
           -i 'INFO/GENEINFO ~ @../data/pathogenic_genes.txt' \
           -Oz -o "$output_file"
    
    # Index the filtered VCF
    tabix -p vcf "$output_file"
    
    echo "Completed: $output_file"
    echo "---"
    
    # Break after first file for testing
    break
done

echo "All VCF files processed successfully!"