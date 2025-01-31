process FASTA2PEPTIDES {
    label 'process_low'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(fasta)

    output:
    tuple val(meta), path("*.tsv"), emit: tsv
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.sample}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII

    """
    fasta2peptides.py \\
        -i $fasta \\
        -o ${prefix} \\
        -minl ${min_length} \\
        -maxl ${max_length}


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.sample}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII

    """
    touch ${prefix}_length_${min_length}.tsv
    touch ${prefix}_length_${max_length}.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
        biopython: \$(python3 -c "import Bio; print(Bio.__version__)")
    END_VERSIONS
    """

}
