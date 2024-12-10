process MHCFLURRY {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcflurry=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.1.4--pyh7e72e81_0' :
        'quay.io/biocontainers/mhcflurry:2.1.4--pyh7e72e81_0' }"

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class == "II") {
        error("MHCflurry prediction of ${meta.sample} is not possible with MHC class II!")
    }
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample

    """
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_mhcflurry.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(mhcflurry-predict --version)
    END_VERSIONS
    """
}
