//and also do the allele name

process PREPARE_PREDICTION_INPUT {
    label 'process_single'
    tag "${meta.sample}"

    conda "bioconda::mhcgnomes=1.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.4--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcgnomes:1.8.4--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.csv"), emit: prepared
    path "versions.yml", emit: versions

    script:
    //TODO handle the thresholds (parse the --tools_thresholds and --use_affinity_thresholds)
    def min_length = (meta.mhc_class == "I") ? params.min_peptide_length_mhc_I : params.min_peptide_length_mhc_II
    def max_length = (meta.mhc_class == "I") ? params.max_peptide_length_mhc_I : params.max_peptide_length_mhc_II
    //tools über params.tools ziehen

    """
    """

    stub:
    """
    touch syfpeithi_input.csv
    touch versions.yml
    """
}
