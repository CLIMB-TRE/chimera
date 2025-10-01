process ALIGNMENT_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/numpy:2.3.3--c8ac2f6d040d90b0'
        : 'community.wave.seqera.io/library/numpy:2.3.3--f39d62792aa0b968'}"

    input:
    tuple val(meta), path(depth_tsv), path(coverage_tsv)
    path database_metadata

    output:
    tuple val(meta), path("*.alignment_report.tsv"), emit: alignment_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_alignment_report.py \\
        --depth_tsv ${depth_tsv} \\
        --coverage_tsv ${coverage_tsv} \\
        --database_metadata ${database_metadata} \\
        > ${prefix}.alignment_report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment_report.tsv
    """
}
