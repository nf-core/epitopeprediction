// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

process DEFINE_SOFTWARE {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'reports', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "bioconda::fred2=2.0.6 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
        if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
            container "https://depot.galaxyproject.org/singularity/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
        } else {
            container "quay.io/biocontainers/mulled-v2-689ae0756dd82c61400782baaa8a7a1c2289930d:a9e10ca22d4cbcabf6b54f0fb5d766ea16bb171e-0"
        }

    output:
       path '*.csv', emit: versions // model_report.txt

    script:
        """
        cat <<-END_VERSIONS > versions.csv
                mhcflurry\t\$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
                mhcnuggets\t\$(echo \$(python -c "import pkg_resources; print('mhcnuggets' + pkg_resources.get_distribution('mhcnuggets').version)" | sed 's/^mhcnuggets//; s/ .*\$//' ))
                fred2\t\$(echo \$(python -c "import pkg_resources; print('fred2' + pkg_resources.get_distribution('Fred2').version)" | sed 's/^fred2//; s/ .*\$//'))
        END_VERSIONS
        """
}


