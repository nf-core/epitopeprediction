// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process CHECK_REQUESTED_MODELS {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'reports', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    } else {
        container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
    }

    input:
        tuple val(alleles), path(input_file)
        val(software_versions)

    output:
        path '*.txt', emit: txt // model_report.txt
        path '*.log', emit: log // model_warnings.log
        val "versions.yml", emit: versions


    script:
        def argument = "$options.args"

        if (argument.contains("peptides") == true) {
            argument += " ${input_file}"
        }

        """
        check_requested_models.py ${argument} \
            --alleles '${alleles.join(';')}' \
            --mhcclass ${params.mhc_class} \
            --versions ${software_versions} > model_warnings.log

        cat <<-END_VERSIONS > versions.yml
        ${getProcessName(task.process)}:
            mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
            mhcnuggets: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)"))
            fred2: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('Fred2').version)"))
        END_VERSIONS
        """

}
