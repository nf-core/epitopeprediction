process MHCFLURRY {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.1.4--pyh7e72e81_1' :
        'quay.io/biocontainers/mhcflurry:2.1.4--pyh7e72e81_1' }"

    // MHCflurry downloads models always to ~/.local/share/mhcflurry
    containerOptions {
        (workflow.containerEngine == 'docker') ? '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME' : ''
    }

    input:
    tuple val(meta), path(csv), path(models)

    output:
    tuple val(meta), path("*.csv"), emit: predicted
    tuple val("${task.process}"), val('mhcflurry'), eval("mhcflurry-predict --version | cut -d' ' -f2"), topic: versions

    script:
    if (meta.mhc_class == "II") {
        error("MHCflurry prediction of ${meta.id} is not possible with MHC class II!")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    export MHCFLURRY_DATA_DIR=$models
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=2.2.0

    mhcflurry-predict \\
        $csv \\
        --out ${prefix}_predicted_mhcflurry.csv \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_mhcflurry.csv
    """
}
