process ALIGNMENT_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_numpy_pysam:410caf1b9aff14b8'
        : 'community.wave.seqera.io/library/pip_numpy_pysam:b0f6802385070dc7'}"

    input:
    tuple val(meta), path(depth_tsv), path(coverage_tsv), path(bam)
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
        --bam ${bam} \\
        > ${prefix}.alignment_report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment_report.tsv
    """
}
