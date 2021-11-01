// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process PEPTIDE_PREDICTION {
    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'', meta:[:], publish_by_meta:[]) }

    // TODO: Change the container
    // fred2:2.0.7-py_0
    // mhcflurry:1.6.1-py_0
    // mhcnuggets:2.3.2-py_0
    // pandas:0.24.2-py27h86efe34_0
    // python:2.7.15-h8e446fc_1011_cpython
    // pyvcf:0.6.8-py27_0
    // snpsift:4.3.1t-hdfd78af_3
    // snpsift=4.3.1t, python=2.7.15, pyvcf=0.6.8, pandas=0.24.2, fred2=2.0.7, mhcflurry=1.6.1, mhcnuggets=2.3.2
    conda (params.enable_conda ? "conda-forge::snpsift:4.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/snpsift:4.2--hdfd78af_5"
    } else {
        container "quay.io/biocontainers/snpsift:4.2--hdfd78af_5"
    }

    input:
        tuple val(meta), path(splitted), val(software_versions) // combine the different paths?
    output:
        tuple val(meta), path("*.tsv"), path("*.json"), path("*.fasta"), emit: predicted

    script:

    """
    epaa.py --identifier ${splitted.baseName} \
                        --alleles "${meta.alleles}" \
                        --versions ${software_versions} \
                        $options.args ${splitted}
    """

}
