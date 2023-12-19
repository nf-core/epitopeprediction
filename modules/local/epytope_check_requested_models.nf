process EPYTOPE_CHECK_REQUESTED_MODELS {
    label 'process_low'

    conda "bioconda::epytope=3.3.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(input_file)
    path(software_versions)

    output:
    path '*.txt', emit: txt // model_report.txt
    path '*.log', emit: log // model_warnings.log
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def argument = task.ext.args

    if (argument.contains("peptides") == true) {
        argument += " ${input_file}"
    }

    if (params.multiqc_title) {
        argument += "--title \"$params.multiqc_title\""
    }

    def prefix = task.ext.suffix ? "${meta.sample}_${task.ext.suffix}" : "${meta.sample}_peptides"
    def min_length = ("${meta.mhc_class}" == "I") ? params.min_peptide_length : params.min_peptide_length_class2
    def max_length = ("${meta.mhc_class}" == "I") ? params.max_peptide_length : params.max_peptide_length_class2

    """
    check_requested_models.py ${argument} \
        --alleles '${meta.alleles}' \
        --max_length ${max_length} \
        --min_length ${min_length} \
        --versions ${software_versions} > model_warnings.log

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mhcflurry: \$(echo \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//') )
        mhcnuggets: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)"))
        epytope: \$(echo \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)"))
    END_VERSIONS
    """
}
