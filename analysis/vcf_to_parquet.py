"""
Export variant-level attributes and per-sample genotype metrics from a VCF into
compressed Parquet tables.

Paths and column requirements are hard-coded so the script can be executed
directly without command-line arguments.
"""

from __future__ import annotations

import os
import math
from pathlib import Path
from typing import Any, Dict, List, Tuple

import polars as pl
from cyvcf2 import VCF


VCF_PATH = Path("output/combined.vcf.bgz")
VARIANT_PARQUET_PATH = Path("output/prevalence/variants.parquet")
SAMPLE_PARQUET_PATH = Path("output/prevalence/samples.parquet")


VARIANT_INFO_SPECS: Dict[str, Tuple[Any, bool, str]] = {
    "CLNSIG": ("", True, "str"),
    "GENEINFO": ("", True, "str"),
    "MC": ("", True, "str"),
    "AF": (0.0, True, "float"),
    "AN": (0, False, "int"),
    "AC": (0, True, "int"),
    "QUALapprox": (0, True, "int"),
    "AS_QUALapprox": (0, True, "int"),
    "SCORE": (0.0, True, "float"),
}


def ensure_parent(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)


def decode_if_bytes(value: Any) -> Any:
    if hasattr(value, "tolist"):
        return decode_if_bytes(value.tolist())
    if isinstance(value, (bytes, bytearray)):
        return value.decode("utf-8")
    if isinstance(value, list):
        return [decode_if_bytes(elem) for elem in value]
    if isinstance(value, tuple):
        return tuple(decode_if_bytes(elem) for elem in value)
    return value


def get_info_value(record, key: str, alt_index: int | None, default: Any) -> Any:
    value = record.INFO.get(key)
    if value is None:
        return default
    value = decode_if_bytes(value)
    if isinstance(value, (list, tuple)):
        if alt_index is not None and len(value) == len(record.ALT):
            return value[alt_index]
        if len(value) == 1:
            return value[0]
        return ",".join("" if elem is None else str(elem) for elem in value)
    return value


def cast_value(value: Any, kind: str, default: Any) -> Any:
    if value in (None, "", ".", []):
        return default
    try:
        if kind == "float":
            return float(value)
        if kind == "int":
            return int(float(value))
        return str(value)
    except (TypeError, ValueError):
        return default


def format_gt(genotype: List[int]) -> str:
    if not genotype:
        return ""
    alleles = genotype[:-1]
    if not alleles:
        return ""
    if all(allele < 0 for allele in alleles):
        return "./."
    phased = bool(genotype[-1])
    tokens = ["." if allele < 0 else str(allele) for allele in alleles]
    separator = "|" if phased else "/"
    return separator.join(tokens)


def format_ad(ad_row: Any) -> List[int]:
    if ad_row is None:
        return []
    ad_row = decode_if_bytes(ad_row)
    if isinstance(ad_row, (int, float)):
        return [int(ad_row)]
    if isinstance(ad_row, str):
        try:
            return [int(float(ad_row))]
        except ValueError:
            return []
    if not isinstance(ad_row, (list, tuple)):
        return []
    formatted: List[int] = []
    for value in ad_row:
        if value in (None, -1, ".", ""):
            formatted.append(0)
            continue
        try:
            formatted.append(int(float(value)))
        except (TypeError, ValueError):
            formatted.append(0)
    return formatted


def scalar_value(array, index: int, default: Any, kind: str = "int") -> Any:
    if array is None:
        return default
    try:
        value = array[index]
    except (IndexError, TypeError):
        return default
    value = decode_if_bytes(value)
    if isinstance(value, (list, tuple)):
        value = value[0] if value else default
    if value in (None, -1, ".", ""):
        return default
    if isinstance(value, float) and math.isnan(value):
        return default
    try:
        if kind == "float":
            return float(value)
        if kind == "str":
            return str(value)
        return int(float(value))
    except (TypeError, ValueError):
        return default


def prepare_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    variant_rows: list[dict[str, Any]] = []
    sample_rows: list[dict[str, Any]] = []
    vcf = VCF(str(VCF_PATH))
    samples = vcf.samples
    try:
        for record in vcf:
            alt_count = len(record.ALT)
            if alt_count == 0:
                continue

            for alt_index, alt in enumerate(record.ALT):
                row: dict[str, Any] = {
                    "chrom": record.CHROM,
                    "pos": record.POS,
                    "ref": record.REF,
                    "alt": "" if alt is None else str(alt),
                }
                for key, (default, per_alt, kind) in VARIANT_INFO_SPECS.items():
                    value = get_info_value(
                        record,
                        key,
                        alt_index if per_alt else None,
                        default,
                    )
                    row[key] = cast_value(value, kind, default)
                variant_rows.append(row)

            genotypes = record.genotypes
            ad_array = record.format("AD")
            gq_array = record.format("GQ")
            ps_array = record.format("PS")
            rgq_array = record.format("RGQ")
            ft_array = record.format("FT")

            for idx, sample in enumerate(samples):
                genotype = list(genotypes[idx]) if genotypes is not None else []
                gt_value = format_gt(genotype)

                ad_value = format_ad(ad_array[idx]) if ad_array is not None else []
                gq_value = scalar_value(gq_array, idx, 0, "int")
                ps_value = scalar_value(ps_array, idx, 0, "int")
                rgq_value = scalar_value(rgq_array, idx, 0, "int")
                ft_value = scalar_value(ft_array, idx, "", "str")

                sample_rows.append(
                    {
                        "chrom": record.CHROM,
                        "pos": record.POS,
                        "ref": record.REF,
                        "alts": ",".join("" if a is None else str(a) for a in record.ALT),
                        "sample": sample,
                        "gt": gt_value,
                        "ad": ad_value,
                        "gq": gq_value,
                        "ps": ps_value,
                        "rgq": rgq_value,
                        "ft": ft_value,
                    }
                )
    finally:
        vcf.close()
    return variant_rows, sample_rows


def write_parquet(rows: list[dict[str, Any]], path: Path) -> None:
    if not rows:
        df = pl.DataFrame()
    else:
        df = pl.DataFrame(rows)
    ensure_parent(path)
    df.write_parquet(path, compression="zstd")


def main() -> None:
    variant_rows, sample_rows = prepare_rows()
    write_parquet(variant_rows, VARIANT_PARQUET_PATH)
    write_parquet(sample_rows, SAMPLE_PARQUET_PATH)


if __name__ == "__main__":
    main()
