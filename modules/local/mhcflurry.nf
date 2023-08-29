process MHCFLURRY {
    label 'process_single'
    tag "${metadata.sample}"

    conda "bioconda::mhcflurry=2.0.6"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcflurry:2.0.6--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcflurry:2.0.6--pyh7cba7a3_0' }"

    input:
    tuple val(metadata), path(peptide_file)

    output:
    tuple val(metadata), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:

    if (metadata.mhc_class == "II") {
        error("MHCflurry prediction of ${metadata.sample} is not possible with MHC class II!")
    }

    def prefix = "${metadata.sample}_${peptide_file.baseName}"
    def min_length = (metadata.mhc_class == "I") ? params.min_peptide_length_mhc_I : params.min_peptide_length_mhc_II
    def max_length = (metadata.mhc_class == "I") ? params.max_peptide_length_mhc_I : params.max_peptide_length_mhc_II

    """
    """

    stub:
    """
    touch ${metadata.sample}_predicted_mhcflurry.tsv
    touch versions.yml
    """
}
