// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SPLIT_PEPTIDES {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'splitted', meta:[:], publish_by_meta:[]) }

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

    script:
        def prefix = options.suffix ? "${peptide.baseName}_${options.suffix}" : "${peptide.baseName}"

        """
        split_peptides.py --input ${peptide} \\
        --output_base "${prefix}" \\
        $options.args
        """
}
