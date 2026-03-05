process PEPTDEEP_LIBRARY {
    tag "$meta.id"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/peptdeep:1.4.1--pyhdfd78af_1' :
        'quay.io/biocontainers/peptdeep:1.4.1--pyhdfd78af_1' }"

    input:
    tuple val(meta), path(peptide_tsv)

    output:
    tuple val(meta), path("*.parquet"), emit: speclib
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    """
    # Writable HOME for peptdeep model downloads (container HOME=/ is read-only)
    mkdir -p nxf_home
    export HOME=\$PWD/nxf_home
    # Limit numpy/pytorch thread backends to allocated CPUs
    export OMP_NUM_THREADS=${task.cpus}
    export OPENBLAS_NUM_THREADS=${task.cpus}
    export MKL_NUM_THREADS=${task.cpus}

    peptdeep_library.py \\
        --input ${peptide_tsv} \\
        --output ${prefix}.parquet \\
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
    touch ${prefix}.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        peptdeep: \$(python -c "import peptdeep; print(peptdeep.__version__)")
        alphabase: \$(python -c "import alphabase; print(alphabase.__version__)")
    END_VERSIONS
    """
}
