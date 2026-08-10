# CLIMB-TRE/chimera: Output

## Introduction

This document describes the output produced by the pipeline.

The directories listed below will be created in the results directory (`--outdir`) after the pipeline has finished, under a subdirectory named after each sample's `id` (as given in the `sample` column of the samplesheet), except for `pipeline_info/` which is shared across the whole run. All paths are relative to the top-level results directory.

## Pipeline overview

The pipeline is built using [Nextflow](https://www.nextflow.io/) and processes data using the following steps, for each sample in the samplesheet:

- [Sylph taxonomic profiling](#sylph-taxonomic-profiling) - Estimate the taxonomic composition of the raw reads against the reference database
- [Read alignment](#read-alignment) - Align reads to the reference database with minimap2 (`ont`) or bwa-mem2 (`illumina`/`illumina.se`)
- [Alignment proportion filtering](#alignment-proportion-filtering) - Optionally remove reads with a low proportion of aligned bases
- [Per-reference depth](#per-reference-depth) - Compute per-base coverage of each reference the sample aligned to
- [Sylph taxonomy report](#sylph-taxonomy-report) - Resolve the sylph profile's reference genomes to NCBI taxIDs and scientific names
- [Alignment report](#alignment-report) - Per-reference alignment statistics and confidence scoring
- [Pipeline information](#pipeline-information) - Report metrics generated during the workflow execution

### Sylph taxonomic profiling

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.tsv`: Raw [sylph](https://github.com/bluenote-1577/sylph) `profile` output, one row per reference genome sylph detected in the sample, with columns including taxonomic/sequence abundance and estimated ANI.

</details>

[sylph](https://github.com/bluenote-1577/sylph) is run against `--sylph_db`, an ANI-based genome/metagenome profiler that estimates which reference genomes are present in a sample and at what abundance, without requiring alignment. This raw profile is later combined with taxonomy information from the [sylph taxonomy report](#sylph-taxonomy-report).

### Read alignment

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.bam`: Coordinate-sorted BAM file of reads aligned against the reference database (`--mm2_index` for `ont` samples via minimap2, `--bwa_index_prefix` for `illumina`/`illumina.se` samples via bwa-mem2).
  - `*.bam.bai`: BAM index.

</details>

ONT reads are aligned with [minimap2](https://github.com/lh3/minimap2) (`-x map-ont --secondary=yes -N 50 --secondary-seq`), which reports up to 50 secondary alignments per read with sequence included. Illumina reads are aligned with [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) (`-a -h 50`); because bwa-mem2 only reports secondary alignments as full records for unpaired reads, alternative loci recorded in each read's `XA` tag are additionally expanded into full secondary alignment records, with sequence and quality copied (reverse-complemented where necessary) from that read's primary alignment, so that Illumina secondary alignments carry sequence in the same way ONT ones do. This expansion is always applied to Illumina alignments and is not user-configurable.

### Alignment proportion filtering

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.bam`: BAM file with reads failing the alignment proportion filter removed (only produced when `--min_alignment_proportion_filter` is `true`; this file replaces the one described in [Read alignment](#read-alignment) above).

</details>

When `--min_alignment_proportion_filter` is enabled, reads (or read pairs, for Illumina data, where either mate passing the threshold is sufficient to retain both) whose proportion of aligned bases is below `--min_alignment_proportion` (default `0.5`) are removed before depth calculation and reporting. This step is skipped entirely by default.

### Per-reference depth

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.tsv`: Per-base depth of coverage for every reference the sample aligned to, as produced by `samtools depth -a`.

</details>

[samtools depth](http://www.htslib.org/doc/samtools-depth.html) is used to compute per-base coverage across every position of every reference sequence the sample's reads aligned to (including zero-coverage positions, via `-a`), which is used to calculate the evenness, breadth and depth metrics in the [alignment report](#alignment-report).

### Sylph taxonomy report

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.sylph_taxonomy_report.tsv`: The [sylph taxonomic profiling](#sylph-taxonomic-profiling) output, joined with the NCBI taxID and scientific name of each reference genome via `--sylph_taxdb`.

</details>

The reference genome accessions in the raw sylph profile are matched against `--sylph_taxdb` (a genome-accession → NCBI-taxID lookup table; see [Building the sylph taxonomy database](usage.md#building-the-sylph-taxonomy-database) in the usage documentation for how to build one) to add taxonomic context to each row of the sylph profile.

### Alignment report

<details markdown="1">
<summary>Output files</summary>

- `<sample_id>/`
  - `*.alignment_report.tsv`: One row per reference genome the sample aligned to, with alignment/coverage statistics and, if `--alignment_scoring_matrix` is used, a `total_score` and `confidence` category.

</details>

For every reference genome a sample's reads aligned to, this report combines:

- Coverage metrics from [per-reference depth](#per-reference-depth): evenness of coverage, mean depth, breadth at 1x/10x, mapped reads/bases.
- Alignment-quality metrics computed directly from the BAM: mean read identity (from the `NM` tag), read duplication rate (reads sharing a start/end position), mean alignment length/proportion, forward-strand proportion, and alignment complexity (a measure of how repetitive the aligned sequence is).
- Reference metadata from `--database_metadata`: taxon ID, human-readable name, accession description and sequence length.

If a scoring matrix (`--alignment_scoring_matrix`, defaulting to [`assets/alignment_scoring_matrix.json`](../assets/alignment_scoring_matrix.json)) is supplied, each metric is scored against a set of configurable bands, summed into a `total_score`, and mapped to a `confidence` category (`fail`, `low_confidence`, `medium_confidence` or `high_confidence` by default) — intended to help distinguish well-supported reference hits from low-confidence or spurious ones. Rows are sorted by descending `total_score`.

### Pipeline information

<details markdown="1">
<summary>Output files</summary>

- `pipeline_info/`
  - Reports generated by Nextflow: `execution_report.html`, `execution_timeline.html`, `execution_trace.txt` and `pipeline_dag.dot`/`pipeline_dag.svg`.
  - Reformatted samplesheet files used as input to the pipeline: `samplesheet.valid.csv`.
  - Parameters used by the pipeline run: `params.json`.
  - Collated software versions: `chimera_software_versions.yml`.

</details>

[Nextflow](https://www.nextflow.io/docs/latest/tracing.html) provides excellent functionality for generating various reports relevant to the running and execution of the pipeline. This will allow you to troubleshoot errors with the running of the pipeline, and also provide you with other information such as launch commands, run times and resource usage.
