process MHCFLURRY {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcflurry=2.1.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.1.4--pyh7e72e81_1' :
        'quay.io/biocontainers/mhcflurry:2.1.4--pyh7e72e81_1' }"

    // MHCflurry downloads models always to ~/.local/share/mhcflurry
    containerOptions = (workflow.containerEngine == 'docker') ? '-u $(id -u) -e "HOME=${HOME}" -v /etc/passwd:/etc/passwd:ro -v /etc/shadow:/etc/shadow:ro -v /etc/group:/etc/group:ro -v $HOME:$HOME' : ''

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.csv"), emit: predicted
    path "versions.yml"           , emit: versions

    script:
    if (meta.mhc_class == "II") {
        error("MHCflurry prediction of ${meta.sample} is not possible with MHC class II!")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.sample}"

    """
    mhcflurry-downloads fetch models_class1_presentation
    mhcflurry-predict \\
        $peptide_file \\
        --out ${prefix}_predicted_mhcflurry.csv \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(mhcflurry-predict --version)
    END_VERSIONS
    """

    stub:
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_mhcflurry.csv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(mhcflurry-predict --version)
    END_VERSIONS
    """
}
