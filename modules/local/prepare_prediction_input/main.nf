process PREPARE_PREDICTION_INPUT {
    label 'process_single'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(tsv)
    path(supported_alleles_json)

    output:
    tuple val(meta), path("*.json"), path("*.{csv,tsv}"), emit: prepared
    path "versions.yml"                                 , emit: versions

    script:
    template "prepare_prediction_input.py"

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mhcflurry_input.csv
    touch ${prefix}_mhcnuggets_input.tsv
    touch ${prefix}_mhcnuggetsii_input.tsv
    touch ${prefix}_netmhcpan_input.tsv
    touch ${prefix}_netmhciipan_input.tsv
    touch versions.yml
    """
}
