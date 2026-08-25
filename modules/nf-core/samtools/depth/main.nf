process SAMTOOLS_DEPTH {
    tag "$meta1.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/samtools:1.22.1--h96c455f_0' :
        'biocontainers/samtools:1.22.1--h96c455f_0' }"

    input:
    tuple val(meta1), path(bam), path(bai)
    tuple val(meta2), path(intervals)

    output:
    tuple val(meta1), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta1.id}"
    def positions = intervals ? "-b ${intervals}" : ""
    """
    if [ -n "${positions}" ]; then
        DEPTH_POSITIONS="${positions}"
    else
        samtools idxstats ${bam} | awk '\$3 > 0 {print \$1"\\t0\\t"\$2}' > covered_regions.bed
        DEPTH_POSITIONS="-b covered_regions.bed"
    fi

    # Note: --threads value represents *additional* CPUs to allocate (total CPUs = 1 + --threads).
    samtools \\
        depth \\
        --threads ${task.cpus-1} \\
        $args \\
        \${DEPTH_POSITIONS} \\
        -o ${prefix}.tsv \\
        $bam

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
