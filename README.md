# CLIMB-TRE/chimera

[![GitHub Actions CI Status](https://github.com/CLIMB-TRE/chimera/actions/workflows/nf-test.yml/badge.svg)](https://github.com/CLIMB-TRE/chimera/actions/workflows/nf-test.yml)
[![GitHub Actions Linting Status](https://github.com/CLIMB-TRE/chimera/actions/workflows/linting.yml/badge.svg)](https://github.com/CLIMB-TRE/chimera/actions/workflows/linting.yml)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/version-%E2%89%A524.10.5-green?style=flat&logo=nextflow&logoColor=white&color=%230DC09D&link=https%3A%2F%2Fnextflow.io)](https://www.nextflow.io/)
[![nf-core template version](https://img.shields.io/badge/nf--core_template-3.3.2-green?style=flat&logo=nfcore&logoColor=white&color=%2324B064&link=https%3A%2F%2Fnf-co.re)](https://github.com/nf-core/tools/releases/tag/3.3.2)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/CLIMB-TRE/chimera)

## Introduction

**CLIMB-TRE/chimera** is a bioinformatics pipeline for profiling Illumina and Nanopore (ONT) sequencing reads against a curated reference database and flagging how confidently each detected reference is actually supported by the read data. It runs [sylph](https://github.com/bluenote-1577/sylph) for fast taxonomic profiling of raw reads, aligns reads to the reference database with [rammap](https://github.com/jwanglab/rammap) (a single minimap2-compatible aligner used for both ONT and Illumina data), and combines per-reference coverage and alignment-quality statistics into a scored confidence report — helping to distinguish genuine hits from cross-mapping, contamination, or other low-confidence alignments.

The pipeline:

1. Profiles each sample's taxonomic composition against a sylph database (`SYLPH_PROFILE`)
2. Aligns reads to a reference database with rammap, using a platform-specific preset (`RAMMAP_ALIGN`)
3. Optionally filters out reads with a low proportion of aligned bases (`FILTER_BAM`)
4. Computes per-base, per-reference depth of coverage (`SAMTOOLS_DEPTH`)
5. Resolves the sylph profile's reference genomes to NCBI taxIDs (`SYLPH_TAXONOMY`)
6. Produces a per-reference alignment report with configurable confidence scoring (`ALIGNMENT_REPORT`)

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

First, prepare a samplesheet with your input data that looks as follows:

`samplesheet.csv`:

```csv
sample,platform,fastq_1,fastq_2
CONTROL_REP1,illumina,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
```

Each row represents a fastq file (`platform` of `ont` or `illumina.se`) or a pair of fastq files (`platform` of `illumina`). See [`docs/usage.md`](docs/usage.md) for the full samplesheet specification and for the reference database parameters (`--rammap_index`, `--sylph_db`, `--sylph_taxdb`, `--database_metadata`) the pipeline requires.

Now, you can run the pipeline using:

```bash
nextflow run CLIMB-TRE/chimera \
   -profile <docker/singularity/.../institute> \
   --input samplesheet.csv \
   --outdir <OUTDIR> \
   --rammap_index /path/to/reference.mmi \
   --sylph_db /path/to/reference.syldb \
   --sylph_taxdb /path/to/sylph_taxdb.tsv \
   --database_metadata /path/to/database_metadata.tsv
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

## Credits

CLIMB-TRE/chimera was originally written by biowilko.

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/main/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
