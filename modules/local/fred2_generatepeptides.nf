process FRED2_GENERATEPEPTIDES {

    conda (params.enable_conda ? "conda-forge::fred2:2.0.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fred2:2.0.7--py_0"
    } else {
        container "quay.io/biocontainers/fred2:2.0.7--py_0"
    }

    input:
        tuple val(meta), path(raw)

    output:
        tuple val(meta), path("*.tsv"), emit: splitted
        path "versions.yml", emit: versions

    script:
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_peptides"

        """
        gen_peptides.py --input ${raw} \\
        --output '${prefix}.tsv' \\
        $options.args

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            fred2: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('Fred2').version)")
            python: \$(python --version 2>&1 | sed 's/Python //g')
        END_VERSIONS
        """
}
