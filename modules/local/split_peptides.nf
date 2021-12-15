process SPLIT_PEPTIDES {

    conda (params.enable_conda ? "conda-forge::python=3.8.3" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/python:3.8.3"
    } else {
        container "quay.io/biocontainers/python:3.8.3"
    }
    // cache false

    input:
        tuple val(meta), path(peptide)

    output:
        tuple val(meta), path("*.tsv"), emit: splitted
        path "versions.yml", emit: versions

    script:
        def prefix = options.suffix ? "${peptide.baseName}_${options.suffix}" : "${peptide.baseName}"

        """
        split_peptides.py --input ${peptide} \\
        --output_base "${prefix}" \\
        $options.args

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            python: \$(python --version | sed 's/Python //g')
        END_VERSIONS
        """
}
