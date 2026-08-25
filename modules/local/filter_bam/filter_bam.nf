process FILTER_BAM {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pysam_jsonschema_numpy:0d4beb5f588f4b4d'
        : 'community.wave.seqera.io/library/pysam_jsonschema_numpy:8ec3d505e678a720'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.filtered.bam"), emit: filtered_bam

    when:
    task.ext.when == null || task.ext.when

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
