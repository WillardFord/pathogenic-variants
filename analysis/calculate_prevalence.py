"""
Compute per-variant carrier prevalence and allele frequency using cyvcf2.

Paths are hard-coded via module constants so the script can be executed
directly without command-line flags.
"""

from __future__ import annotations

import csv
from pathlib import Path
from typing import Iterator, Sequence

from cyvcf2 import VCF


# Edit these constants to point at the desired inputs/outputs.
VCF_PATH = Path("output/combined_filtered.vcf.bgz")
OUTPUT_PATH = Path("output/prevalence/by_variant.tsv")


Row = Sequence[object]


def prevalence_rows(vcf_path: Path) -> Iterator[Row]:
    """Yield prevalence rows for each alternate allele in the VCF/BCF."""
    vcf = VCF(str(vcf_path))
    total_samples = len(vcf.samples)
    try:
        for record in vcf:
            alt_count = len(record.ALT)
            if alt_count == 0:
                continue

            carriers = [0] * alt_count
            alt_copies = [0] * alt_count
            called_samples = 0
            total_alleles = 0

            for genotype in record.genotypes:
                if not genotype:
                    continue
                alleles = genotype[:-1]  # last element encodes phasing
                valid = [allele for allele in alleles if allele >= 0]
                if not valid:
                    continue
                called_samples += 1
                total_alleles += len(valid)
                alt_indices = set()
                for allele in valid:
                    if allele > 0:
                        alt_copies[allele - 1] += 1
                        alt_indices.add(allele)
                for allele_index in alt_indices:
                    if 1 <= allele_index <= alt_count:
                        carriers[allele_index - 1] += 1

            for offset, alt in enumerate(record.ALT):
                carrier_count = carriers[offset]
                total_samples = (
                    total_samples
                )
                # Called samples should be identical to total samples
                called_samples = (
                    called_samples
                )
                allele_count = alt_copies[offset]
                # TODO fix below
                yield [
                    record.CHROM,
                    record.POS,
                    record.REF,
                    alt,
                    carrier_count,
                    total_samples,
                    called_samples,
                    allele_count,
                ]
    finally:
        vcf.close()


def write_rows(rows: Iterator[Row], output_path: Path) -> None:
    header = [
        "chrom",
        "pos",
        "ref",
        "alt",
        "carrier_count",
        "total_samples",
        "called_samples",
        "allele_count",
    ]
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(header)
        for row in rows:
            writer.writerow(row)


def main() -> None:
    write_rows(prevalence_rows(VCF_PATH), OUTPUT_PATH)


if __name__ == "__main__":
    main()

"""
For every sample we need
GT
AD
GQ
PS
RGQ
FT

For every variant we need:
CLNSIG
GENEINFO
MC
AF
AN
AC
QUALapprox
AS_QUALapprox
SCORE
"""