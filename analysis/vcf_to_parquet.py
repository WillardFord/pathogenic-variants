"""
Export variant-level attributes and per-sample genotype metrics from VCF shards.

When annotated shard files exist under
``output/filtered_vcfs_with_overlapping_ranges/``, each file is converted in
parallel into separate Parquet pieces beneath
``output/prevalence/parquet/{variants,samples}``. If no shards are found, the
script falls back to the combined VCF at ``output/combined.vcf.bgz`` and writes
single Parquet tables for variants and samples. Paths and column requirements
remain hard-coded so the script can be executed without command-line arguments.
"""

from __future__ import annotations

import math
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Tuple

import polars as pl
import pyarrow.parquet as pq
from cyvcf2 import VCF


VCF_PATH = Path("output/combined.vcf.bgz")
VARIANT_PARQUET_PATH = Path("output/prevalence/variants.parquet")
SAMPLE_PARQUET_PATH = Path("output/prevalence/samples.parquet")
ANNOTATED_VCF_DIR = Path("output/filtered_vcfs_with_overlapping_ranges")
VARIANT_DATASET_DIR = Path("output/prevalence/parquet/variants")
SAMPLE_DATASET_DIR = Path("output/prevalence/parquet/samples")


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


def positive_int_from_env(name: str, default: int) -> int:
    value = os.environ.get(name)
    if value is None:
        return default
    try:
        parsed = int(value)
    except ValueError as exc:
        raise ValueError(f"{name} must be a positive integer (got {value!r})") from exc
    if parsed <= 0:
        raise ValueError(f"{name} must be greater than zero (got {parsed})")
    return parsed


def bool_from_env(name: str, default: bool) -> bool:
    value = os.environ.get(name)
    if value is None:
        return default
    normalized = value.strip().lower()
    if normalized in {"1", "true", "yes", "on"}:
        return True
    if normalized in {"0", "false", "no", "off"}:
        return False
    raise ValueError(f"{name} must be a truthy/falsey string (got {value!r})")


@dataclass
class ConversionConfig:
    vcf_threads: int
    records_per_chunk: int
    variant_rows_per_chunk: int
    sample_rows_per_chunk: int
    parquet_threads: int


def build_config_from_env() -> ConversionConfig:
    config = ConversionConfig(
        vcf_threads=positive_int_from_env("VCF_TO_PARQUET_VCF_THREADS", 16),
        records_per_chunk=positive_int_from_env("VCF_TO_PARQUET_RECORDS_PER_CHUNK", 100),
        variant_rows_per_chunk=positive_int_from_env("VCF_TO_PARQUET_VARIANT_ROWS_PER_CHUNK", 2000),
        sample_rows_per_chunk=positive_int_from_env("VCF_TO_PARQUET_SAMPLE_ROWS_PER_CHUNK", 50_000),
        parquet_threads=positive_int_from_env("VCF_TO_PARQUET_PARQUET_THREADS", 16),
    )
    return config


def configure_thread_env(config: ConversionConfig) -> None:
    os.environ.setdefault("PYARROW_NUM_THREADS", str(config.parquet_threads))
    os.environ.setdefault("POLARS_MAX_THREADS", str(config.parquet_threads))


def strip_known_suffixes(name: str) -> str:
    lowercase = name.lower()
    for suffix in (".vcf.bgz", ".vcf.gz", ".bgz", ".gz", ".vcf"):
        if lowercase.endswith(suffix):
            return name[: -len(suffix)]
    return Path(name).stem


def dataset_output_paths(vcf_path: Path) -> tuple[Path, Path]:
    base = strip_known_suffixes(vcf_path.name)
    variant_path = VARIANT_DATASET_DIR / f"{base}.parquet"
    sample_path = SAMPLE_DATASET_DIR / f"{base}.parquet"
    return variant_path, sample_path


def find_annotated_vcfs(directory: Path) -> List[Path]:
    if not directory.is_dir():
        return []
    candidates = set(directory.glob("*_annotated.vcf.bgz"))
    candidates.update(directory.glob("*_annotated.vcf.gz"))
    return sorted(candidates)


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


def empty_variant_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": pl.Series("chrom", [], dtype=pl.Utf8),
            "pos": pl.Series("pos", [], dtype=pl.Int64),
            "ref": pl.Series("ref", [], dtype=pl.Utf8),
            "alt": pl.Series("alt", [], dtype=pl.Utf8),
            "CLNSIG": pl.Series("CLNSIG", [], dtype=pl.Utf8),
            "GENEINFO": pl.Series("GENEINFO", [], dtype=pl.Utf8),
            "MC": pl.Series("MC", [], dtype=pl.Utf8),
            "AF": pl.Series("AF", [], dtype=pl.Float64),
            "AN": pl.Series("AN", [], dtype=pl.Int64),
            "AC": pl.Series("AC", [], dtype=pl.Int64),
            "QUALapprox": pl.Series("QUALapprox", [], dtype=pl.Int64),
            "AS_QUALapprox": pl.Series("AS_QUALapprox", [], dtype=pl.Int64),
            "SCORE": pl.Series("SCORE", [], dtype=pl.Float64),
        }
    )


def empty_sample_df() -> pl.DataFrame:
    return pl.DataFrame(
        {
            "chrom": pl.Series("chrom", [], dtype=pl.Utf8),
            "pos": pl.Series("pos", [], dtype=pl.Int64),
            "ref": pl.Series("ref", [], dtype=pl.Utf8),
            "alts": pl.Series("alts", [], dtype=pl.Utf8),
            "sample": pl.Series("sample", [], dtype=pl.Utf8),
            "gt": pl.Series("gt", [], dtype=pl.Utf8),
            "ad": pl.Series("ad", [], dtype=pl.List(pl.Int64)),
            "gq": pl.Series("gq", [], dtype=pl.Int64),
            "ps": pl.Series("ps", [], dtype=pl.Int64),
            "rgq": pl.Series("rgq", [], dtype=pl.Int64),
            "ft": pl.Series("ft", [], dtype=pl.Utf8),
        }
    )


class ParquetStreamWriter:
    def __init__(self, path: Path, empty_frame_factory) -> None:
        self._path = path
        self._writer: pq.ParquetWriter | None = None
        self._empty_frame_factory = empty_frame_factory
        ensure_parent(path)

    def write(self, rows: list[dict[str, Any]]) -> None:
        if not rows:
            return
        df = pl.DataFrame(rows)
        if df.is_empty():
            return
        table = df.to_arrow()
        if self._writer is None:
            self._writer = pq.ParquetWriter(self._path, table.schema, compression="zstd")
        self._writer.write_table(table)

    def close(self) -> None:
        if self._writer is not None:
            self._writer.close()
            self._writer = None
        else:
            empty_df = self._empty_frame_factory()
            empty_df.write_parquet(self._path, compression="zstd")


def iter_vcf_batches(vcf_path: Path, config: ConversionConfig) -> Iterable[tuple[list[dict[str, Any]], list[dict[str, Any]]]]:
    vcf = VCF(str(vcf_path))
    vcf.set_threads(config.vcf_threads)
    samples = vcf.samples
    variant_batch: list[dict[str, Any]] = []
    sample_batch: list[dict[str, Any]] = []
    record_counter = 0

    try:
        for record in vcf:
            alt_count = len(record.ALT)
            if alt_count == 0:
                continue

            record_counter += 1

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
                variant_batch.append(row)

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

                sample_batch.append(
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

            should_flush = False
            if config.records_per_chunk and record_counter % config.records_per_chunk == 0:
                should_flush = True
            if len(variant_batch) >= config.variant_rows_per_chunk:
                should_flush = True
            if len(sample_batch) >= config.sample_rows_per_chunk:
                should_flush = True

            if should_flush:
                yield variant_batch, sample_batch
                variant_batch, sample_batch = [], []

        if variant_batch or sample_batch:
            yield variant_batch, sample_batch
    finally:
        vcf.close()


def convert_vcf_to_parquet(
    vcf_path: Path,
    variant_path: Path,
    sample_path: Path,
    config: ConversionConfig | None = None,
    force: bool = False,
) -> None:
    if config is None:
        config = build_config_from_env()
    configure_thread_env(config)

    if not force and variant_path.exists() and sample_path.exists():
        print(f"[SKIP] {vcf_path.name}: outputs already exist.", flush=True)
        return

    print(f"[CONVERT] {vcf_path.name} -> {variant_path.name}", flush=True)

    variant_writer = ParquetStreamWriter(variant_path, empty_variant_df)
    sample_writer = ParquetStreamWriter(sample_path, empty_sample_df)

    try:
        for variant_rows, sample_rows in iter_vcf_batches(vcf_path, config):
            variant_writer.write(variant_rows)
            sample_writer.write(sample_rows)
    finally:
        variant_writer.close()
        sample_writer.close()

    print(f"[DONE] {vcf_path.name}", flush=True)


def _convert_worker(
    vcf_path: str,
    variant_path: str,
    sample_path: str,
    config: ConversionConfig,
    force: bool,
) -> str:
    convert_vcf_to_parquet(
        Path(vcf_path),
        Path(variant_path),
        Path(sample_path),
        config=config,
        force=force,
    )
    return vcf_path


def main() -> None:
    config = build_config_from_env()
    force = bool_from_env("VCF_TO_PARQUET_FORCE", False)
    default_procs = max(os.cpu_count() or 1, 1)
    max_processes = positive_int_from_env("VCF_TO_PARQUET_MAX_PROCESSES", default_procs)

    annotated_vcfs = find_annotated_vcfs(ANNOTATED_VCF_DIR)

    if annotated_vcfs:
        VARIANT_DATASET_DIR.mkdir(parents=True, exist_ok=True)
        SAMPLE_DATASET_DIR.mkdir(parents=True, exist_ok=True)

        pending = []
        for vcf_path in annotated_vcfs:
            variant_path, sample_path = dataset_output_paths(vcf_path)
            if not force and variant_path.exists() and sample_path.exists():
                print(f"[SKIP] {vcf_path.name}: outputs already exist.")
                continue
            pending.append((vcf_path, variant_path, sample_path))

        if not pending:
            print("All annotated VCFs already converted; nothing to do.")
            return

        worker_count = min(len(pending), max_processes)
        print(f"Converting {len(pending)} VCF shard(s) with {worker_count} process(es)...")

        with ProcessPoolExecutor(max_workers=worker_count) as pool:
            futures = [
                pool.submit(
                    _convert_worker,
                    str(vcf_path),
                    str(variant_path),
                    str(sample_path),
                    config,
                    force,
                )
                for vcf_path, variant_path, sample_path in pending
            ]
            for future in as_completed(futures):
                future.result()
        return

    convert_vcf_to_parquet(
        VCF_PATH,
        VARIANT_PARQUET_PATH,
        SAMPLE_PARQUET_PATH,
        config=config,
        force=force,
    )


if __name__ == "__main__":
    main()
