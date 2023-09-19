process MHCNUGGETS {
    label 'process_low'
    tag "${meta.sample}"

    conda "bioconda::mhcnuggets=2.4.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcnuggets:2.4.1--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcnuggets:2.4.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    def prefix = "${meta.sample}_${peptide_file.baseName}"
    """
    """

    stub:
    """
    touch ${meta.sample}_predicted_mhcnuggets.tsv
    touch versions.yml
    """
}
