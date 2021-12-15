process PEPTIDE_PREDICTION {
    label 'process_low'

    conda (params.enable_conda ? "bioconda::snpsift=4.3.1t bioconda::python=2.7.15 bioconda::pyvcf=0.6.8 conda-forge::pandas=0.24.2 bioconda::fred2=2.0.7 bioconda::mhcflurry=1.4.3 bioconda::mhcnuggets=2.3.2" : null)
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0' :
        'quay.io/biocontainers/mulled-v2-c3f301504f7fa2e7bf81c3783de19a9990ea3001:12b1b9f040fd92a80629d58f8a558dde4820eb15-0' }"

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
        fred2: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('Fred2').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyvcf').version)")
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}
