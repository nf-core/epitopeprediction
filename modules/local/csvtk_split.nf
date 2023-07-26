process CSVTK_SPLIT {
    label 'process_low'

    conda "conda-forge::sed=4.7 bioconda::csvtk=0.23.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(raw)

    output:
    tuple val(meta), path("*.tsv"), emit: splitted
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    sed -i.bak '/^##/d' ${raw}
    csvtk split ${raw} -t -C '&' -f '#chr'

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """

}
