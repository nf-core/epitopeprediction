process SUMMARIZE_RESULTS {
    label 'process_low'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.tsv") , emit: tsv
    tuple val(meta), path("*.json"), emit: json
    tuple val("${task.process}"), val('python'), eval("python --version | cut -d' ' -f2"), topic: versions, emit: versions_python

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_results.py \\
        --input . \\
        --prefix ${prefix} \\
        --peptide_col_name ${params.peptide_col_name} \\
        $args

    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    touch ${prefix}_mqc.json
    """
}
