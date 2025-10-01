#!/bin/bash

curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz -o data/clinvar.vcf.gz
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz.tbi -o data/clinvar.vcf.gz.tbi
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_papu.vcf.gz -o data/clinvar_papu.vcf.gz
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_papu.vcf.gz.tbi -o data/clinvar_papu.vcf.gz.tbi

# The vcfs we use need to have 'chr' prefix on chromosomes to match All of Us data
echo "Adding 'chr' prefix to ClinVar chromosomes..."
bcftools annotate --rename-chrs ../data/chr_mapping.txt -Oz -o ../data/clinvar_with_chr.vcf.gz ../data/clinvar.vcf.gz
tabix -p vcf ../data/clinvar_with_chr.vcf.gz

echo "Adding 'chr' prefix to ClinVar PaPu chromosomes..."
bcftools annotate --rename-chrs ../data/chr_mapping.txt -Oz -o ../data/clinvar_papu_with_chr.vcf.gz ../data/clinvar_papu.vcf.gz
tabix -p vcf ../data/clinvar_papu_with_chr.vcf.gz


# And genocde gtf
curl https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz -o data/gencode.v49.primary_assembly.annotation.gtf.gz