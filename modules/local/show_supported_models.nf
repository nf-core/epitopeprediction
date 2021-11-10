// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process SHOW_SUPPORTED_MODELS {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'supported_models', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(meta), path(raw), val(software_versions)

    output:
        path '*.txt', emit: txt // model_report.txt

    script:

        """
        check_supported_models.py --versions ${software_versions}
        """

}
