"""
Extract unique GENEINFO annotations from the VCF and write them to disk.
"""

from __future__ import annotations

from pathlib import Path

from cyvcf2 import VCF


VCF_PATH = Path("output/combined.vcf.gz")
OUTPUT_PATH = Path("output/prevalence/geneinfo_tags.txt")


def decode_if_bytes(value):
    if hasattr(value, "tolist"):
        return decode_if_bytes(value.tolist())
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8")
    if isinstance(value, list):
        return [decode_if_bytes(elem) for elem in value]
    if isinstance(value, tuple):
        return tuple(decode_if_bytes(elem) for elem in value)
    return value


def main() -> None:
    tags = set()
    vcf = VCF(str(VCF_PATH))
    try:
        for record in vcf:
            value = record.INFO.get("GENEINFO")
            if value is None:
                continue
            value = decode_if_bytes(value)
            if isinstance(value, (list, tuple)):
                for elem in value:
                    if elem:
                        tags.add(str(elem))
            else:
                if value:
                    tags.add(str(value))
    finally:
        vcf.close()

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    with OUTPUT_PATH.open("w") as handle:
        for tag in sorted(tags):
            handle.write(f"{tag}\n")


if __name__ == "__main__":
    main()
