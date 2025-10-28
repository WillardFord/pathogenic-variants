"""
Count unique individuals carrying at least one alternate allele in the VCF.

Paths are hard-coded so the script can be executed directly.
"""

from __future__ import annotations
from collections import defaultdict
from pathlib import Path
from typing import Any

from cyvcf2 import VCF


VCF_PATH = Path("output/gene_filtered/filtered_pathogenic_vcfs/combined.vcf.gz")
OUTPUT_PATH = Path("output/prevalence/carrier_count.txt")
PATHOGENICITY_REQUIREMENTS_PATH = Path("data/ACMG SF Annotations Download - ACMG SF v3.3 List.tsv")

# Pre-load gene pathogenicity requirements from the file into a dictionary for efficient lookup.
def load_pathogenicity_requirements(path: Path) -> dict[str, str]:
    requirements = {}
    with path.open() as handle:
        for line in handle:
            if line.startswith("Gene"):
                continue
            fields = line.rstrip('\n').split('\t')
            # ACMG gene symbol is field[0], requirement is field[7]
            requirements[fields[0].strip().lower()] = fields[7].strip()
    return requirements

GENE_PATHOGENICITY_REQUIREMENTS = load_pathogenicity_requirements(PATHOGENICITY_REQUIREMENTS_PATH)

def get_pathogenicity_requirements(geneinfo: str) -> str:
    """Return the pathogenicity requirements for a gene."""
    # Try to find gene requirement by substring matching
    for gene, req in GENE_PATHOGENICITY_REQUIREMENTS.items():
        if gene in geneinfo.lower():
            return req
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
    carriers: dict[str, list[str]] = defaultdict(list)
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
                    carriers[samples[idx]].append(record.INFO.get("GENEINFO")+"_het")
                elif pathogenic == True:
                    carriers[samples[idx]].append(record.INFO.get("GENEINFO"))

    finally:
        vcf.close()

    for sample, genes in carriers.items():
        cleaned_genes = []
        gene_set = set(genes)
        for gene in genes:
            if gene.endswith("_het") and genes.count(gene) < 2:
                continue
            else:
                cleaned_genes.append(gene)
        if len(cleaned_genes) == 0:
            del carriers[sample]
        carriers[sample] = cleaned_genes

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    count = len(carriers)

    # TODO change carriers
    with OUTPUT_PATH.open("w") as handle:
        handle.write(f"carrier_count\t{count}\n")
        for sample, genes in carriers.items():
            handle.write(f"{sample}\t{'\t'.join(genes)}\n")


if __name__ == "__main__":
    main()
