# CLIMB-TRE/chimera: Usage

> _Documentation of pipeline parameters is generated automatically from the pipeline schema and can no longer be found in markdown files._

## Introduction

CLIMB-TRE/chimera takes Illumina or Nanopore (ONT) sequencing reads, profiles their taxonomic composition with [sylph](https://github.com/bluenote-1577/sylph), and aligns them against a reference database of candidate taxa (with [minimap2](https://github.com/lh3/minimap2) for ONT reads or [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2) for Illumina reads) to produce a per-reference alignment quality report, scored against a configurable confidence matrix. This is intended to help distinguish genuine, well-supported alignments to a reference taxon from spurious or low-confidence hits (e.g. cross-mapping between related taxa, contamination, or chimeric reads) flagged by the taxonomic profiler.

## Samplesheet input

You will need to create a samplesheet with information about the samples you would like to analyse before running the pipeline. Use the `--input` parameter to specify its location. It has to be a comma-separated file with 4 columns (3 required), and a header row as shown in the examples below.

```bash
--input '[path to samplesheet file]'
```

### Multiple runs of the same sample

The `sample` identifiers have to be the same when you have re-sequenced the same sample more than once e.g. to increase sequencing depth. The pipeline will concatenate the raw reads before performing any downstream analysis. Below is an example for the same sample sequenced across 3 lanes:

```csv title="samplesheet.csv"
sample,platform,fastq_1,fastq_2
CONTROL_REP1,illumina,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP1,illumina,AEG588A1_S1_L003_R1_001.fastq.gz,AEG588A1_S1_L003_R2_001.fastq.gz
CONTROL_REP1,illumina,AEG588A1_S1_L004_R1_001.fastq.gz,AEG588A1_S1_L004_R2_001.fastq.gz
```

### Full samplesheet

The pipeline uses the `platform` column, not file count, to decide how a sample is treated: `illumina` (paired-end), `illumina.se` (single-end Illumina) or `ont` (Nanopore, always single-end). The samplesheet can have as many columns as you desire, however, there is a strict requirement for the first 4 columns to match those defined in the table below.

A final samplesheet file consisting of Illumina paired-end, Illumina single-end and ONT data may look something like the one below.

```csv title="samplesheet.csv"
sample,platform,fastq_1,fastq_2
CONTROL_REP1,illumina,AEG588A1_S1_L002_R1_001.fastq.gz,AEG588A1_S1_L002_R2_001.fastq.gz
CONTROL_REP2,illumina,AEG588A2_S2_L002_R1_001.fastq.gz,AEG588A2_S2_L002_R2_001.fastq.gz
TREATMENT_REP1,illumina.se,AEG588A4_S4_L003_R1_001.fastq.gz,
TREATMENT_REP2,ont,AEG588A5_S5_L003_R1_001.fastq.gz,
```

| Column     | Description                                                                                                                                                                              |
| ---------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `sample`   | Custom sample name. This entry will be identical for multiple sequencing libraries/runs from the same sample. Spaces in sample names are automatically converted to underscores (`_`). |
| `platform` | Sequencing platform used for this sample. Must be one of `illumina` (paired-end), `illumina.se` (single-end Illumina), or `ont` (Nanopore).                                             |
| `fastq_1`  | Full path to a FastQ file. File has to be gzipped and have the extension ".fastq.gz" or ".fq.gz".                                                                                        |
| `fastq_2`  | Full path to the second read of an `illumina` paired-end FastQ pair. Required when `platform` is `illumina`; must be left empty for `illumina.se` and `ont` samples.                    |

An [example samplesheet](../assets/samplesheet.csv) has been provided with the pipeline.

## Reference database inputs

Unlike pipelines that build their own reference index, CLIMB-TRE/chimera expects a pre-built reference database and its associated indexes/metadata to be supplied via parameters, so the (potentially large) reference-building step can be done once and reused across many runs:

| Parameter                        | Description                                                                                                                                                              |
| --------------------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `--mm2_index`                     | Path to a pre-built minimap2 index (`.mmi`) of the reference database, used for `ont` samples.                                                                          |
| `--bwa_index_prefix`               | Path prefix to a pre-built bwa-mem2 index of the reference database (e.g. `ref.fa` if the index files are `ref.fa.0123`, `ref.fa.amb`, etc.), used for Illumina samples. |
| `--sylph_db`                       | Path to a pre-sketched sylph database (`.syldb`) of the reference genomes, used for taxonomic profiling.                                                                |
| `--sylph_taxdb`                    | Path to a sylph taxonomy lookup TSV mapping each reference genome accession to an NCBI taxID and scientific name. See [Building the sylph taxonomy database](#building-the-sylph-taxonomy-database) below for how to build this. |
| `--database_metadata`              | Path to a TSV with one row per reference accession, containing at minimum `unique_accession`, `taxon_id`, `human_readable`, `accession_description` and `sequence_length` columns, used to annotate the alignment report. |
| `--alignment_scoring_matrix`       | Path to the JSON file defining the banded scoring matrix used to score each reference's alignment statistics (defaults to the bundled [`assets/alignment_scoring_matrix.json`](../assets/alignment_scoring_matrix.json)). |
| `--alignment_scoring_json_schema`  | Path to the JSON schema used to validate `--alignment_scoring_matrix` (defaults to the bundled [`assets/alignment_scoring_matrix_schema.json`](../assets/alignment_scoring_matrix_schema.json)). |

All reference genomes referenced by `--bwa_index_prefix`/`--mm2_index` should be the same set as those sketched into `--sylph_db`, keyed by the same accessions used in `--sylph_taxdb` and `--database_metadata`, so that results can be joined across the sylph and alignment arms of the pipeline.

## Running the pipeline

The typical command for running the pipeline is as follows:

```bash
nextflow run CLIMB-TRE/chimera \
    --input ./samplesheet.csv \
    --outdir ./results \
    --mm2_index /path/to/reference.mmi \
    --bwa_index_prefix /path/to/reference.fa \
    --sylph_db /path/to/reference.syldb \
    --sylph_taxdb /path/to/sylph_taxdb.tsv \
    --database_metadata /path/to/database_metadata.tsv \
    -profile docker
```

This will launch the pipeline with the `docker` configuration profile. See below for more information about profiles.

Note that the pipeline will create the following files in your working directory:

```bash
work                # Directory containing the nextflow working files
<OUTDIR>            # Finished results in specified location (defined with --outdir)
.nextflow_log       # Log file from Nextflow
# Other nextflow hidden files, eg. history of pipeline runs and old logs.
```

If you wish to repeatedly use the same parameters for multiple runs, rather than specifying each flag in the command, you can specify these in a params file.

Pipeline settings can be provided in a `yaml` or `json` file via `-params-file <file>`.

> [!WARNING]
> Do not use `-c <file>` to specify parameters as this will result in errors. Custom config files specified with `-c` must only be used for [tuning process resource specifications](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources), other infrastructural tweaks (such as output directories), or module arguments (args).

The above pipeline run specified with a params file in yaml format:

```bash
nextflow run CLIMB-TRE/chimera -profile docker -params-file params.yaml
```

with:

```yaml title="params.yaml"
input: './samplesheet.csv'
outdir: './results/'
<...>
```

You can also generate such `YAML`/`JSON` files via [nf-core/launch](https://nf-co.re/launch).

### Updating the pipeline

When you run the above command, Nextflow automatically pulls the pipeline code from GitHub and stores it as a cached version. When running the pipeline after this, it will always use the cached version if available - even if the pipeline has been updated since. To make sure that you're running the latest version of the pipeline, make sure that you regularly update the cached version of the pipeline:

```bash
nextflow pull CLIMB-TRE/chimera
```

### Reproducibility

It is a good idea to specify the pipeline version when running the pipeline on your data. This ensures that a specific version of the pipeline code and software are used when you run your pipeline. If you keep using the same tag, you'll be running the same version of the pipeline, even if there have been changes to the code since.

First, go to the [CLIMB-TRE/chimera releases page](https://github.com/CLIMB-TRE/chimera/releases) and find the latest pipeline version - numeric only (eg. `1.3.1`). Then specify this when running the pipeline with `-r` (one hyphen) - eg. `-r 1.3.1`. Of course, you can switch to another version by changing the number after the `-r` flag.

This version number will be logged in reports when you run the pipeline, so that you'll know what you used when you look back in the future.

To further assist in reproducibility, you can use share and reuse [parameter files](#running-the-pipeline) to repeat pipeline runs with the same settings without having to write out a command with every single parameter.

> [!TIP]
> If you wish to share such profile (such as upload as supplementary material for academic publications), make sure to NOT include cluster specific paths to files, nor institutional specific profiles.

## Alignment proportion filtering

By default the pipeline does **not** filter reads by alignment proportion. Set `--min_alignment_proportion_filter` to `true` to enable filtering, and `--min_alignment_proportion` (default `0.5`) to set the minimum proportion of a read that must be aligned for it to be kept. When filtering is disabled, `--min_alignment_proportion` is ignored.

## Illumina secondary alignments

Illumina alignments are produced with bwa-mem2 run with `-a -h 50`, and multi-mapping loci recorded in `XA` tags are expanded into full secondary alignment records (SEQ/QUAL copied, with reverse-complementing as needed, from the read's primary alignment) so that Illumina secondary alignments carry sequence in the same way ONT secondary alignments already do via minimap2's `--secondary-seq`.

## Building the sylph taxonomy database

`params.sylph_taxdb` is not produced by the pipeline itself — it is built offline, ahead of time, with `bin/sylph_ncbi_taxonomy_rationaliser.py`. This script takes the raw sylph reference taxonomy (`contig_fname`, `taxon_string` columns) and resolves each entry to an NCBI taxID by joining on GTDB assembly accession against a GTDB metadata TSV (e.g. `bac120_metadata_r232.tsv.gz` from the [GTDB downloads page](https://gtdb.ecogenomic.org/downloads)), falling back to a literal (non-fuzzy) match of the GTDB organism name/lineage against NCBI taxonomy names and synonyms.

```bash
bin/sylph_ncbi_taxonomy_rationaliser.py \
    /path/to/ncbi_taxdump_dir \
    /path/to/raw_sylph_taxonomy.tsv \
    --gtdb_metadata /path/to/bac120_metadata_r232.tsv.gz \
    --report /path/to/unresolved_taxa_report.tsv \
    > sylph_taxdb.tsv
```

The taxdump directory must contain `names.dmp` and `nodes.dmp` (and, if available, `merged.dmp`/`delnodes.dmp` so merged/deleted taxIDs are handled correctly). `--gtdb_metadata` can be given more than once (e.g. to add `ar53_metadata_r232.tsv.gz` for archaea).

This tool is interactive: whenever a taxon cannot be resolved via a direct accession/taxID match, it prompts on the terminal (reading from `/dev/tty`, so stdin/stdout stay free for the input/output taxonomy files) with any literal name matches it found, and asks for confirmation before accepting one. Every unresolved case — whether confirmed or skipped — is written incrementally to the `--report` TSV as the run progresses, so a long run can be safely interrupted without losing the report so far.

## Core Nextflow arguments

> [!NOTE]
> These options are part of Nextflow and use a _single_ hyphen (pipeline parameters use a double-hyphen)

### `-profile`

Use this parameter to choose a configuration profile. Profiles can give configuration presets for different compute environments.

Several generic profiles are bundled with the pipeline which instruct the pipeline to use software packaged using different methods (Docker, Singularity, Podman, Shifter, Charliecloud, Apptainer, Conda) - see below.

> [!IMPORTANT]
> We highly recommend the use of Docker or Singularity containers for full pipeline reproducibility, however when this is not possible, Conda is also supported.

The pipeline also dynamically loads configurations from [https://github.com/nf-core/configs](https://github.com/nf-core/configs) when it runs, making multiple config profiles for various institutional clusters available at run time. For more information and to check if your system is supported, please see the [nf-core/configs documentation](https://github.com/nf-core/configs#documentation).

Note that multiple profiles can be loaded, for example: `-profile test,docker` - the order of arguments is important!
They are loaded in sequence, so later profiles can overwrite earlier profiles.

If `-profile` is not specified, the pipeline will run locally and expect all software to be installed and available on the `PATH`. This is _not_ recommended, since it can lead to different results on different machines dependent on the computer environment.

- `test`
  - A profile with a complete configuration for automated testing
  - Includes links to test data so needs no other parameters
- `docker`
  - A generic configuration profile to be used with [Docker](https://docker.com/)
- `singularity`
  - A generic configuration profile to be used with [Singularity](https://sylabs.io/docs/)
- `podman`
  - A generic configuration profile to be used with [Podman](https://podman.io/)
- `shifter`
  - A generic configuration profile to be used with [Shifter](https://nersc.gitlab.io/development/shifter/how-to-use/)
- `charliecloud`
  - A generic configuration profile to be used with [Charliecloud](https://hpc.github.io/charliecloud/)
- `apptainer`
  - A generic configuration profile to be used with [Apptainer](https://apptainer.org/)
- `wave`
  - A generic configuration profile to enable [Wave](https://seqera.io/wave/) containers. Use together with one of the above (requires Nextflow ` 24.03.0-edge` or later).
- `conda`
  - A generic configuration profile to be used with [Conda](https://conda.io/docs/). Please only use Conda as a last resort i.e. when it's not possible to run the pipeline with Docker, Singularity, Podman, Shifter, Charliecloud, or Apptainer.

### `-resume`

Specify this when restarting a pipeline. Nextflow will use cached results from any pipeline steps where the inputs are the same, continuing from where it got to previously. For input to be considered the same, not only the names must be identical but the files' contents as well. For more info about this parameter, see [this blog post](https://www.nextflow.io/blog/2019/demystifying-nextflow-resume.html).

You can also supply a run name to resume a specific run: `-resume [run-name]`. Use the `nextflow log` command to show previous run names.

### `-c`

Specify the path to a specific config file (this is a core Nextflow command). See the [nf-core website documentation](https://nf-co.re/usage/configuration) for more information.

## Custom configuration

### Resource requests

Whilst the default requirements set within the pipeline will hopefully work for most people and with most input data, you may find that you want to customise the compute resources that the pipeline requests. Each step in the pipeline has a default set of requirements for number of CPUs, memory and time. For most of the pipeline steps, if the job exits with any of the error codes specified [here](https://github.com/nf-core/rnaseq/blob/4c27ef5610c87db00c3c5a3eed10b1d161abf575/conf/base.config#L18) it will automatically be resubmitted with higher resources request (2 x original, then 3 x original). If it still fails after the third attempt then the pipeline execution is stopped.

To change the resource requests, please see the [max resources](https://nf-co.re/docs/usage/configuration#max-resources) and [tuning workflow resources](https://nf-co.re/docs/usage/configuration#tuning-workflow-resources) section of the nf-core website.

### Custom Containers

In some cases, you may wish to change the container or conda environment used by a pipeline steps for a particular tool. By default, nf-core pipelines use containers and software from the [biocontainers](https://biocontainers.pro/) or [bioconda](https://bioconda.github.io/) projects. However, in some cases the pipeline specified version maybe out of date.

To use a different container from the default container or conda environment specified in a pipeline, please see the [updating tool versions](https://nf-co.re/docs/usage/configuration#updating-tool-versions) section of the nf-core website.

### Custom Tool Arguments

A pipeline might not always support every possible argument or option of a particular tool used in pipeline. Fortunately, nf-core pipelines provide some freedom to users to insert additional parameters that the pipeline does not include by default.

To learn how to provide additional arguments to a particular tool of the pipeline, please see the [customising tool arguments](https://nf-co.re/docs/usage/configuration#customising-tool-arguments) section of the nf-core website.

### nf-core/configs

In most cases, you will only need to create a custom config as a one-off but if you and others within your organisation are likely to be running nf-core pipelines regularly and need to use the same settings regularly it may be a good idea to request that your custom config file is uploaded to the `nf-core/configs` git repository. Before you do this please can you test that the config file works with your pipeline of choice using the `-c` parameter. You can then create a pull request to the `nf-core/configs` repository with the addition of your config file, associated documentation file (see examples in [`nf-core/configs/docs`](https://github.com/nf-core/configs/tree/master/docs)), and amending [`nfcore_custom.config`](https://github.com/nf-core/configs/blob/master/nfcore_custom.config) to include your custom profile.

See the main [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) for more information about creating your own configuration files.

If you have any questions or issues please send us a message on [Slack](https://nf-co.re/join/slack) on the [`#configs` channel](https://nfcore.slack.com/channels/configs).

## Running in the background

Nextflow handles job submissions and supervises the running jobs. The Nextflow process must run until the pipeline is finished.

The Nextflow `-bg` flag launches Nextflow in the background, detached from your terminal so that the workflow does not stop if you log out of your session. The logs are saved to a file.

Alternatively, you can use `screen` / `tmux` or similar tool to create a detached session which you can log back into at a later time.
Some HPC setups also allow you to run nextflow within a cluster job submitted your job scheduler (from where it submits more jobs).

## Nextflow memory requirements

In some cases, the Nextflow Java virtual machines can start to request a large amount of memory.
We recommend adding the following line to your environment to limit this (typically in `~/.bashrc` or `~./bash_profile`):

```bash
NXF_OPTS='-Xms1g -Xmx4g'
```
