process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcgnomes=1.8.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(prediction_files)

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    path "versions.yml", emit: versions


    script:
    //TODO handle the thresholds (parse the --tools_thresholds and --use_affinity_thresholds)
    template "merge_predictions.py"

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch merged_prediction.tsv
    touch versions.yml
    """
}
