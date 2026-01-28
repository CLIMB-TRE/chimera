process SYLPH_TAXONOMY {
    tag "${meta.id}"
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'oras://community.wave.seqera.io/library/pysam_jsonschema_numpy:0d4beb5f588f4b4d'
        : 'community.wave.seqera.io/library/pysam_jsonschema_numpy:8ec3d505e678a720'}"

    input:
    tuple val(meta), path(sylph_profile)
    path sylph_taxonomy_lookup

    output:
    tuple val(meta), path("*.sylph_taxonomy_report.tsv"), emit: sylph_taxonomy_report

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_sylph_taxonomy_report.py \\
        --sylph_report ${sylph_profile} \\
        --sylph_taxonomy ${sylph_taxonomy_lookup} \\
        > ${prefix}.sylph_taxonomy_report.tsv
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.sylph_taxonomy_report.tsv
    """
}
