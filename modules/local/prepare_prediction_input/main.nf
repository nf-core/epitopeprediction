process PREPARE_PREDICTION_INPUT {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.6--pyh7cba7a3_0' :
        'biocontainers/mhcgnomes:1.8.6--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(tsv)
    path(supported_alleles_json)

    output:
    tuple val(meta), path("*_allele_input.json"), path("*_input.{csv,tsv}", arity: '1..*'), emit: prepared // arity: a single file must still arrive as a list
    path "versions.yml", topic: versions

    script:
    template "prepare_prediction_input.py"

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def inputs = params.tools.tokenize(',').collect { tool ->
        def ext = tool == 'mhcflurry' ? 'csv' : 'tsv'
        [tool: tool, filename: "${prefix}_${tool}_input.${ext}"]
    }
    def manifest = inputs.collect { i -> """{"tool": "${i.tool}", "alleles": "HLA-A*01:01", "chunk_id": "", "alleles_input": "HLA-A*01:01", "filename": "${i.filename}"}""" }.join(',')
    """
    touch ${inputs.collect { i -> i.filename }.join(' ')}
    echo '[${manifest}]' > ${prefix}_allele_input.json
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //')
    END_VERSIONS
    """
}
