process PEPTIDE_PREDICTION {
    label 'process_low'

    conda (params.enable_conda ? "conda-forge::gawk=5.1.0 bioconda::epytope=3.0.0" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.0.0--pyh5e36f6f_0' :
        'quay.io/biocontainers/epytope:3.0.0--pyh5e36f6f_0' }"

    input:
    tuple val(meta), path(splitted), path(software_versions)

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.tsv"), emit: predicted
    tuple val(meta), path("*.fasta"), emit: fasta optional true
    path "versions.yml", emit: versions

    script:
    // Additions to the argument command need to go to the beginning.
    // Argument list needs to end with --peptides or --somatic_mutation
    def argument = task.ext.args

    if (params.proteome) {
        argument = "--proteome ${params.proteome} " + argument
    }

    if (params.wild_type) {
        argument = "--wild_type " + argument
    }

    if (params.fasta_output) {
        argument = "--fasta_output " + argument
    }

    if (params.tool_thresholds) {
        argument = "--tool_thresholds ${params.tool_thresholds} " + argument
    }

    """
    epaa.py --identifier ${splitted.baseName} \
        --alleles '${meta.alleles}' \
        --versions ${software_versions} \
        ${argument} ${splitted}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyvcf').version)")
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}
