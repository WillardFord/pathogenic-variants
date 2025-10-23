"""
Count unique individuals carrying at least one alternate allele in the VCF.

Paths are hard-coded so the script can be executed directly.
"""

from __future__ import annotations

from pathlib import Path

from cyvcf2 import VCF
from typing import Any

VCF_PATH = Path("output/gene_filtered/filtered_pathogenic_vcfs/combined.vcf.gz")
OUTPUT_PATH = Path("output/prevalence/carrier_count.txt")
PATHOGENICITY_REQUIREMENTS_PATH = Path("data/ACMG SF Annotations Download - ACMG SF v3.3 List.tsv")

def get_pathogenicity_requirements(geneinfo: str) -> str:
    """Return the pathogenicity requirements for a gene."""
    with PATHOGENICITY_REQUIREMENTS_PATH.open() as handle:
        for line in handle:
            # Skip the header
            if line.startswith("Gene"): continue
            fields = line.rstrip('\n').split('\t')
            if fields[0].lower() not in geneinfo.lower(): continue
            return fields[7]
    raise ValueError(f"No pathogenicity requirements found for {geneinfo}")

def pathogenic_variant(alleles: list[int], geneinfo: str, row: dict[str, Any]) -> bool:
    gene_pathogenicity_requirements = get_pathogenicity_requirements(geneinfo)
    match gene_pathogenicity_requirements:
        case "All hemi or homozygous P and LP or 2 het. P/LP variants": # TODO this should include 2 different variants
            count = sum(1 for allele in alleles if allele > 0)
            if count == 1 and len(alleles) == 1:
                return True
            elif count == 1 and len(alleles) > 1:
                return "Het"
            elif count == 2: # homozygous
                return True
            else:
                return False
        case "All hemi, het, homozygous P and LP":
            return any(allele > 0 for allele in alleles)
        case "All P and LP":
            return any(allele > 0 for allele in alleles)
        case "P and LP (2 variants)": 
            return sum(1 for allele in alleles if allele > 0) == 2
        case "P and LP (truncating variants only)":
            truncating_variants = ["SO:0001575|splice_donor_variant", "SO:0001574|splice_acceptor_variant", "SO:0001587|stop_gained", "SO:0001589|frameshift_variant"]
            if any(variant in row["mc"] for variant in truncating_variants):
                return any(allele > 0 for allele in alleles)
            else:
                return False
        case "p.C282Y homozygotes only":
            if row["CHROM"] == "chr6" and row["POS"] == 26092913:
                return any(allele > 0 for allele in alleles)
            else:
                return False
        case _:
            return False

def main() -> None:
    carriers = set()
    vcf = VCF(str(VCF_PATH))
    samples = vcf.samples
    try:
        for record in vcf:
            genotypes = record.genotypes
            if genotypes is None:
                continue
            for idx, genotype in enumerate(genotypes):

                if not genotype:
                    continue
                alleles = genotype[:-1]
                if not alleles:
                    continue
                pathogenic = pathogenic_variant(alleles, record.INFO.get("GENEINFO"), record)
                if pathogenic == "Het":
                    pass
                elif pathogenic == True:
                    carriers.add(samples[idx])

    finally:
        vcf.close()

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    count = len(carriers)
    with OUTPUT_PATH.open("w") as handle:
        handle.write(f"carrier_count\t{count}\n")
        for sample in sorted(carriers):
            handle.write(f"{sample}\n")


if __name__ == "__main__":
    main()
