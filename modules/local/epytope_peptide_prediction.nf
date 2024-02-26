process EPYTOPE_PEPTIDE_PREDICTION {
    label 'process_low'

    conda "conda-forge::coreutils=9.1 conda-forge::tcsh=6.20.00 bioconda::epytope=3.1.0 conda-forge::gawk=5.1.0 conda-forge::perl=5.32.1"
    container 'ghcr.io/jonasscheid/epitopeprediction-2:0.3.0'

    input:
    tuple val(meta), path(splitted), path(software_versions)
    val netmhc_paths

    output:
    tuple val(meta), path("*.json"), emit: json
    tuple val(meta), path("*.tsv"), emit: predicted, optional: true
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

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

    if (params.use_affinity_thresholds) {
        argument = "--use_affinity_thresholds " + argument
    }

    def netmhc_paths_string = netmhc_paths.join(",")
    def tools_split = params.tools.split(',')
    def class1_tools = tools_split.findAll { ! it.matches('.*(?i)(class-2|ii).*') }
    def class2_tools = tools_split.findAll { it.matches('.*(?i)(syf|class-2|ii).*') }

    if (((meta.mhc_class == "I") & class1_tools.empty) | ((meta.mhc_class == "II") & class2_tools.empty)) {
        exit 1, "No tools specified for mhc class ${meta.mhc_class}"
    }

    def min_length = (meta.mhc_class == "I") ? params.min_peptide_length : params.min_peptide_length_class2
    def max_length = (meta.mhc_class == "I") ? params.max_peptide_length : params.max_peptide_length_class2

    def tools_to_use = ((meta.mhc_class == "I") | (meta.mhc_class == "H-2")) ? class1_tools.join(',') : class2_tools.join(',')

    """
    # create folder for MHCflurry downloads to avoid permission problems when running pipeline with docker profile and mhcflurry selected
    mkdir -p mhcflurry-data
    export MHCFLURRY_DATA_DIR=./mhcflurry-data
    # specify MHCflurry release for which to download models, need to be updated here as well when MHCflurry will be updated
    export MHCFLURRY_DOWNLOADS_CURRENT_RELEASE=1.4.0
    # Add non-free software to the PATH
    shopt -s nullglob
    IFS=',' read -r -a netmhc_paths_string <<< \"$netmhc_paths_string\"
    for p in "\${netmhc_paths_string[@]}"; do
            export PATH="\$(realpath -s "\$p"):\$PATH";
        done
    shopt -u nullglob

    epaa.py --identifier ${splitted.baseName} \
        --alleles '${meta.alleles}' \
        --tools '${tools_to_use}' \
        --max_length ${max_length} \
        --min_length ${min_length} \
        --versions ${software_versions} \
        ${argument} ${splitted}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('PyVCF3').version)")
        mhcflurry: \$(mhcflurry-predict --version 2>&1 | sed 's/^mhcflurry //; s/ .*\$//')
        mhcnuggets: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcnuggets').version)")
    END_VERSIONS
    """
}
