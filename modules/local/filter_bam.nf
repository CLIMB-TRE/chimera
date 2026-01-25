process FILTER_BAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pip_numpy_pysam:410caf1b9aff14b8'
        : 'community.wave.seqera.io/library/pip_numpy_pysam:b0f6802385070dc7'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.filtered.bam"), emit: filtered_bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    bam_filter.py \\
        --min_alignment_proportion ${params.min_alignment_proportion} \\
        ${bam} \\
        ${prefix}.filtered.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.filtered.bam
    """
}
