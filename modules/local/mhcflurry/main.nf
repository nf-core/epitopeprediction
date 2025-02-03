process MHCFLURRY {
    label 'process_single'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.1.4--pyh7e72e81_1' :
        'quay.io/biocontainers/mhcflurry:2.1.4--pyh7e72e81_1' }"

    // MHCflurry downloads models always to ~/.local/share/mhcflurry
    containerOptions = (workflow.containerEngine == 'docker') ? '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME' : ''

    input:
    tuple val(meta), path(csv)

    output:
    tuple val(meta), path("*.csv"), emit: predicted
    path "versions.yml"           , emit: versions

script:
    if (meta.mhc_class == "II") {
        error("MHCflurry prediction of ${meta.sample} is not possible with MHC class II!")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    # Create MHCflurry data directory to avoid permission issues
    mkdir -p mhcflurry-data
    export MHCFLURRY_DATA_DIR=./mhcflurry-data
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=2.2.0

    # Check if models are already available
    if ! mhcflurry-downloads info | grep -qE '\\bYES\\b'; then
        mhcflurry-downloads fetch models_class1_presentation
    fi

    mhcflurry-predict \\
        $csv \\
        --out ${prefix}_predicted_mhcflurry.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(mhcflurry-predict --version)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_mhcflurry.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(mhcflurry-predict --version)
    END_VERSIONS
    """
}
