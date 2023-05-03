process GET_PREDICTION_VERSIONS {
    label 'process_low'

    conda "bioconda::epytope=3.1.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.1.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/epytope:3.1.0--pyh5e36f6f_0' }"

    input:
    val external_tool_versions

    output:
    path "versions.csv", emit: versions // this versions.csv is needed by a downstream process as well. TODO: fix to use versions.yml

    script:
    def external_tools = external_tool_versions.join(",")

    """
    cat <<-END_VERSIONS > versions.csv
    mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
    mhcnuggets: \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' )
    epytope: \$(python -c "import pkg_resources; print('epytope' + pkg_resources.get_distribution('epytope').version)" | sed 's/^epytope//; s/ .*\$//')
    END_VERSIONS

    IFS=',' read -r -a external_tools <<< \"$external_tools\"
    if ! [ -z "${external_tool_versions}" ]; then
        for TOOL in "\${external_tools[@]}"; do
            echo "\$TOOL" >> versions.csv
        done
    fi
    """
}


