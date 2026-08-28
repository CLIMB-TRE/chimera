process RAMMAP_ALIGN {
    tag "$meta.id"
    label 'process_high'

    // No conda environment yet - rammap has no bioconda recipe, container only.
    container 'quay.io/biowilko/rammap_samtools:latest'

    input:
    tuple val(meta), path(reads)
    tuple val(meta2), path(reference)
    val bam_index_extension

    output:
    tuple val(meta), path("*.bam")                       , emit: bam
    tuple val(meta), path("*.bam.${bam_index_extension}"), optional: true, emit: index
    path "versions.yml"                                  , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args  = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    def bam_index = bam_index_extension ? "${prefix}.bam##idx##${prefix}.bam.${bam_index_extension} --write-index" : "${prefix}.bam"

    """
    rammap \\
        $args \\
        -a \\
        -t $task.cpus \\
        $reference \\
        $reads \\
        | samtools sort -@ ${task.cpus-1} -o $bam_index $args2 -

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>&1 | sed 's/^rammap //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.bam
    touch ${prefix}.bam.${bam_index_extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        rammap: \$(rammap --version 2>&1 | sed 's/^rammap //')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
