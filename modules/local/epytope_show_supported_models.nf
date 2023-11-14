process EPYTOPE_SHOW_SUPPORTED_MODELS {
    label 'process_low'

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(raw), path(software_versions)

    output:
    path '*.txt', emit: txt // model_report.txt
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    check_supported_models.py --versions ${software_versions}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
        mhcnuggets: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)"))
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
    END_VERSIONS
    """

}
