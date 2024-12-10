process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcgnomes=1.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.4--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcgnomes:1.8.4--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(prediction_files)

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    path "versions.yml", emit: versions

    script:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample

    """
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch merged_prediction.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        VERSION_PLACEHOLDER
    END_VERSIONS
    """

}
