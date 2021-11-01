// Import generic module functions
include { initOptions; saveFiles; getSoftwareName; getProcessName } from './functions'

params.options = [:]
options        = initOptions(params.options)

def VERSION = "2.0.7"

process GEN_PEPTIDES {

    publishDir "${params.outdir}",
        mode: params.publish_dir_mode,
        saveAs: { filename -> saveFiles(filename:filename, options:params.options, publish_dir:'', meta:[:], publish_by_meta:[]) }

    conda (params.enable_conda ? "conda-forge::fred2:2.0.7" : null)
    if (workflow.containerEngine == 'singularity' && !params.singularity_pull_docker_container) {
        container "https://depot.galaxyproject.org/singularity/fred2:2.0.7--py_0"
    } else {
        container "quay.io/biocontainers/fred2:2.0.7--py_0"
    }

    input:
        tuple val(meta), path(raw)

    output:
        tuple val(meta), path("*.tsv"), emit: splitted

    script:
        def prefix = options.suffix ? "${meta.sample}_${options.suffix}" : "${meta.sample}_peptides"

        """
        gen_peptides.py --input ${raw} \\
        --output '${prefix}.tsv' \\
        $options.args

        cat <<-END_VERSIONS > versions.yml
            ${getProcessName(task.process)}:
                fred2: \$(echo \$(python -c "import pkg_resources; print 'fred2 ' + pkg_resources.get_distribution('Fred2').version" | sed 's/^fred2 //; s/ .*\$//') )
            END_VERSIONS
        """
}
