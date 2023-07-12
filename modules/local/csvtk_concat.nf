process CSVTK_CONCAT {
    label 'process_low'

    conda "bioconda::csvtk=0.23.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/csvtk:0.23.0--h9ee0642_0' :
        'biocontainers/csvtk:0.23.0--h9ee0642_0' }"

    input:
    tuple val(meta), path(predicted)

    output:
    tuple val(meta), path("*prediction_result.tsv"), emit: predicted
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    csvtk concat -t $predicted > ${meta.sample}_prediction_result_unsorted.tmp
    csvtk sort -k chr:n,length:n ${meta.sample}_prediction_result_unsorted.tmp -t --out-file ${meta.sample}_prediction_result.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        csvtk: \$(echo \$( csvtk version | sed -e "s/csvtk v//g" ))
    END_VERSIONS
    """
}
