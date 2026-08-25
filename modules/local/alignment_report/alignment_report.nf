process ALIGNMENT_REPORT {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pysam_jsonschema_numpy:0d4beb5f588f4b4d'
        : 'community.wave.seqera.io/library/pysam_jsonschema_numpy:8ec3d505e678a720'}"

    input:
    tuple val(meta), path(depth_tsv), path(bam)
    path database_metadata
    path scoring_matrix
    path json_schema

    output:
    tuple val(meta), path("*.alignment_report.tsv"), emit: alignment_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_alignment_report.py \\
        --depth_tsv ${depth_tsv} \\
        --database_metadata ${database_metadata} \\
        --scoring_matrix ${scoring_matrix} \\
        --json_schema ${json_schema} \\
        ${bam} \\
        > ${prefix}.alignment_report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.alignment_report.tsv
    """
}
