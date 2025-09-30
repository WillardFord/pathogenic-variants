#!/bin/bash

curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz -o data/clinvar.vcf.gz
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar.vcf.gz.tbi -o data/clinvar.vcf.gz.tbi
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_papu.vcf.gz -o data/clinvar_papu.vcf.gz
curl https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/weekly/clinvar_papu.vcf.gz.tbi -o data/clinvar_papu.vcf.gz.tbi