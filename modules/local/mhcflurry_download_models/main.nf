process MHCFLURRY_DOWNLOAD_MODELS {
    label 'process_single'

    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.1.4--pyh7e72e81_1' :
        'quay.io/biocontainers/mhcflurry:2.1.4--pyh7e72e81_1' }"

    output:
    path "mhcflurry-data", emit: models

    script:
    """
    export MHCFLURRY_DATA_DIR=./mhcflurry-data
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=2.2.0
    mhcflurry-downloads fetch models_class1_presentation
    """

    stub:
    """
    mkdir -p mhcflurry-data
    """
}
