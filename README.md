
## Automated ClinVar Filtering Script

The repository ships with `src/filter_vcf_by_gene.sh`, a resumable pipeline that batches ClinVar shards, screens their Picard `interval_list` files against the gene target list, and then annotates only the overlapping VCFs. Key behavior:

- Uses `gsutil ls` plus `gsutil -m cp` to stage interval lists in batches (`BATCH_SIZE`, default 300) with parallel threads (`GSUTIL_THREADS`).
- Builds a combined interval list per batch, tagging each record with its shard ID, then runs `IntervalListTools --ACTION OVERLAPS` against `data/pathogenic_genes_1000000bp.interval_list` to find overlapping shards quickly.
- Honors `START_INDEX` for both interval screening and VCF processing so reruns can resume midstream.
- When a BED-filtered BCF has no variants, emits a header-only VCF and placeholder `.tbi`, increments `skipped`, and moves on without attempting ClinVar annotation.

Run from the repo root with environment variables set (at minimum `GOOGLE_PROJECT`; optionally override `START_INDEX`, `BATCH_SIZE`, or `GSUTIL_THREADS`).
