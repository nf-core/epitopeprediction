process PVACSEQ_INSTALL_VEP_PLUGIN {
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pvactools:7.0.1--pyhdfd78af_0'
        : 'biocontainers/pvactools:7.0.1--pyhdfd78af_0'}"

    input:
    val trigger

    output:
    path "vep_plugins/*.pm", emit: plugins
    tuple val("${task.process}"), val('pvactools'), eval("pip show pvactools | grep '^Version:' | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Taken from the pinned pvactools package, so the plugins always match the container using them.
    """
    mkdir -p vep_plugins
    pvacseq install_vep_plugin vep_plugins
    """

    stub:
    """
    mkdir -p vep_plugins
    touch vep_plugins/Wildtype.pm vep_plugins/Frameshift.pm
    """
}
