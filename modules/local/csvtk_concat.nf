process CSVTK_CONCAT {

    conda (params.enable_conda ? "bioconda::csvtk=0.23.0" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0"
    } else {
        container "quay.io/biocontainers/csvtk:0.23.0--h9ee0642_0"
    }

    input:
        tuple val(meta), path(predicted)

    output:
        tuple val(meta), path("*.tsv"), emit: predicted
        path "versions.yml", emit: versions

    script:
    """
    csvtk concat -t $predicted > ${meta.sample}_prediction_result.tsv

    cat <<-END_VERSIONS > versions.yml
    ${getProcessName(task.process)}:
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
