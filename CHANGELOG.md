# CLIMB-TRE/chimera: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v0.4.0 - 2026-08-10

### `Added`

- `min_alignment_proportion_filter` parameter to make BAM alignment-proportion filtering optional (defaults to `false`)
- Expansion of bwa-mem2 `XA` tags into secondary alignment records with SEQ/QUAL filled in, so Illumina secondary alignments carry sequence like ONT ones already do (`FILL_SECONDARY_SEQ` module, `bin/fill_secondary_seq.py`)
- GTDB metadata accession join in `bin/sylph_ncbi_taxonomy_rationaliser.py`, with interactive confirmation of literal name matches for unresolved taxa and a TSV report of all unresolved cases
- `skip_sylph` parameter to skip sylph profiling and taxonomy reporting entirely (`--sylph_db`/`--sylph_taxdb` not required when enabled)
- Support for an optional `segment` column in `--database_metadata` (e.g. for multi-segment viral genomes); included in the alignment report when at least one reference has a non-blank value

### `Fixed`

- Fixed `read_pair_generator()`/`read_checks()` in `bin/bam_filter.py` dropping or incorrectly pairing secondary alignments
- Fixed `generate_alignment_report.py` crashing on zero-length alignments (e.g. seq-less secondary records)
- Fixed `bin/fill_secondary_seq.py` writing corrupt BAM records (CIGAR/SEQ length mismatch) when a bwa-mem2 `XA` tag's CIGAR didn't match the primary alignment's read length
- Fixed `bin/fill_secondary_seq.py` running out of memory (OOM) on large BAMs by name-collating the input with `samtools collate` and processing one query name's records at a time, instead of buffering every primary and secondary alignment in the file in memory
- Removed fuzzy (Levenshtein/BK-tree) matching from `bin/sylph_ncbi_taxonomy_rationaliser.py` in favour of exact accession/taxID/name matching

### `Changed`

- Renamed alignment report confidence labels in `assets/alignment_scoring_matrix.json` from `low_confidence`/`medium_confidence`/`high_confidence` to `low`/`medium`/`high`

### `Dependencies`

### `Deprecated`
