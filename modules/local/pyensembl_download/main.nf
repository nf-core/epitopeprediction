process PYENSEMBL_DOWNLOAD {
    label 'process_single'
    tag "${genome_reference}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:4.0.1--pyhdfd78af_1' :
        'biocontainers/epytope:4.0.1--pyhdfd78af_1' }"

    input:
    val(genome_reference)

    output:
    path("pyensembl_cache"), emit: cache
    path "versions.yml"    , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    """
    epaa.py \\
        --download_cache \\
        --genome_reference ${genome_reference} \\
        --pyensembl_cache_dir pyensembl_cache \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "from importlib.metadata import version; print(version('epytope'))")
        pyensembl: \$(python -c "from importlib.metadata import version; print(version('pyensembl'))")
    END_VERSIONS
    """

    stub:
    """
    mkdir -p pyensembl_cache

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "from importlib.metadata import version; print(version('epytope'))")
        pyensembl: \$(python -c "from importlib.metadata import version; print(version('pyensembl'))")
    END_VERSIONS
    """
}
