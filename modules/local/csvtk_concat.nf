// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CSVTK_CONCAT {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'predictions', meta:[:], publish_by_meta:[]) }

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
