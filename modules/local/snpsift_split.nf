// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SNPSIFT_SPLIT {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'.', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::snpsift:4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snpsift:4.2--hdfd78af_5"
    } else {
        container "quay.io/biocontainers/snpsift:4.2--hdfd78af_5"
    }

    input:
        tuple val(meta), path(input_file)

    output:
        tuple val(meta), path("*.vcf"), emit: splitted
        path "versions.yml", emit: versions

    script:
    """
    SnpSift split ${input_file}

    cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            snpsift: \$(echo \$(snpsift -version 2>&1 | sed -n 3p | cut -d\$' ' -f3))
    END_VERSIONS
    """
}
