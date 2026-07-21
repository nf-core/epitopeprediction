process VARIANT_FASTA2PEPTIDES {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(annotated_fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII
    def wild_type  = params.wild_type ? '--wild-type' : ''
    """
    variant_fasta2peptides.py \\
        --in-fasta ${annotated_fasta} \\
        --output-prefix ${prefix} \\
        --min-length ${min_length} \\
        --max-length ${max_length} \\
        --peptide-col-name ${params.peptide_col_name} \\
        ${wild_type}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII
    """
    touch ${prefix}_length_${min_length}.tsv
    touch ${prefix}_length_${max_length}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
