process SYFPEITHI {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::epytope=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.1.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/epytope:3.1.0--pyh5e36f6f_0' }"

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
    touch ${meta.sample}_predicted_syfpeithi.tsv
    touch versions.yml
    """

}
