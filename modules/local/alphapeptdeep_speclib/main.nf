process ALPHAPEPTDEEP_SPECLIB {
    tag "$meta.id"
    label 'process_medium'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/pip_peptdeep:1.4.1--e7b1a24eed0dbe66' :
        'community.wave.seqera.io/library/pip_peptdeep:1.4.1--4f3d825268498c74' }"

    input:
    tuple val(meta), path(peptide_tsv)

    output:
    tuple val(meta), path("*.speclib.tsv"), emit: speclib
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    generate_speclib.py \\
        --input ${peptide_tsv} \\
        --output ${prefix}.speclib.tsv \\
        --peptide_col_name ${params.peptide_col_name} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peptdeep: \$(python -c "import peptdeep; print(peptdeep.__version__)")
        alphabase: \$(python -c "import alphabase; print(alphabase.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.speclib.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peptdeep: \$(python -c "import peptdeep; print(peptdeep.__version__)")
        alphabase: \$(python -c "import alphabase; print(alphabase.__version__)")
    END_VERSIONS
    """
}
