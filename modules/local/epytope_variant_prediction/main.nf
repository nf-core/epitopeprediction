process EPYTOPE_VARIANT_PREDICTION {
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/epytope:3.3.1--pyh7cba7a3_0' :
        'biocontainers/epytope:3.3.1--pyh7cba7a3_0' }"

    input:
    tuple val(meta), path(vcf)
    path(biomart_dump)

    output:
    tuple val(meta), path("*.tsv")  , emit: tsv
    tuple val(meta), path("*.fasta"), emit: fasta, optional: true
    path "versions.yml"             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args
    def prefix = task.ext.prefix ?: "${meta.id}"
    def biomart = params.biomart_dump_path ? "--biomart_dump ${biomart_dump}" : ""
    def min_length = (meta.mhc_class == "I") ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = (meta.mhc_class == "I") ? params.max_peptide_length_classI : params.max_peptide_length_classII
    def flanking_region_size = params.fasta_peptide_flanking_region_size

    """
    epaa.py \
        -i ${vcf} \
        -p ${prefix} \
        ${biomart} \
        --max_length ${max_length} \
        --min_length ${min_length} \
        --flanking_region_size ${flanking_region_size} \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('PyVCF3').version)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.tsv
    touch ${prefix}.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        epytope: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('epytope').version)")
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        pyvcf: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('PyVCF3').version)")
    END_VERSIONS
    """
}
