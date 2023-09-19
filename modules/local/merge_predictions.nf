process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcgnomes=1.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.4--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcgnomes:1.8.4--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(prediction_files), path(metadata_file)

    output:
    tuple val(meta), path("*.tsv"), emit: merged
    path "versions.yml", emit: versions

    script:
    def output = prediction_files.first().baseName.split("_").dropRight(2).join("_")

    """
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python \$(python --version | sed 's/Python //g')
        mhcgnomes \$(python -c "from mhcgnomes import version; print(version.__version__)"   )
    END_VERSIONS
    """

    stub:
    """
    touch merged_prediction.tsv
    touch versions.yml
    """

}
