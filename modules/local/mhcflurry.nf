process MHCFLURRY {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcflurry=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.0.6--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcflurry:2.0.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:

    if (meta.mhc_class == "II") {
        error("MHCflurry prediction of ${meta.sample} is not possible with MHC class II!")
    }

    def prefix = "${meta.sample}_${peptide_file.baseName}"

    """
    """

    stub:
    """
    touch ${meta.sample}_predicted_mhcflurry.tsv
    touch versions.yml
    """
}
