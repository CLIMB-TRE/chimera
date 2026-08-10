process FILL_SECONDARY_SEQ {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pysam_jsonschema_numpy:0d4beb5f588f4b4d'
        : 'community.wave.seqera.io/library/pysam_jsonschema_numpy:8ec3d505e678a720'}"

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("*.secondary_filled.bam"), emit: bam

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    fill_secondary_seq.py \\
        ${bam} \\
        ${prefix}.secondary_filled.bam
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.secondary_filled.bam
    """
}
