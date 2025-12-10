
bcftools view -h ../data/clinvar.vcf.gz > ../data/acmg_variants.vcf
gene_pat=$(paste -sd'|' ../data/ar_genes.txt); bcftools view -H ../data/clinvar.vcf.gz | rg -v 'Conflicting' | rg 'Pathogenic|Likely_pathogenic' | rg -F -f ../data/pathogenic_genes.txt | rg -v 'KREMEN1' | rg -P "^(?!.*(?:$gene_pat)).*$|^(?=.*(?:$gene_pat)).*(?:splice_donor_variant|splice_acceptor_variant|stop_gained|frameshift_variant)" | rg -P "^(?!.*HFE:).*$|^(?=.*HFE:)(?=.*26092913).*" >> ../data/acmg_variants.vcf

bcftools annotate --rename-chrs ../data/chr_mapping.txt -Oz -o ../data/acmg_with_chr.vcf.gz ../data/acmg_variants.vcf
bgzip -c ../data/acmg_with_chr.vcf > ../data/acmg_with_chr.vcf.gz
tabix -p vcf ../data/acmg_with_chr.vcf.gz