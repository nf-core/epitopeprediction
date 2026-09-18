process FASTA2PEPTIDES {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(fasta), path(variants_tsv)
    path proteome_reference

    output:
    tuple val(meta), path("*.tsv")           , emit: tsv
    tuple val(meta), path("*.annotated.fasta"), emit: annotated_fasta, optional: true
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions, emit: versions_python
    tuple val("${task.process}"), val('biopython'), eval('python3 -c "import Bio; print(Bio.__version__)"'), topic: versions, emit: versions_biopython

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII
    def variant    = variants_tsv ? "--variants-tsv ${variants_tsv} --annotated-fasta ${prefix}.annotated.fasta" : ''
    def wild_type  = variants_tsv && params.wild_type ? '--wild-type' : ''
    def proteome   = variants_tsv && proteome_reference ? "--proteome-reference ${proteome_reference}" : ''
    """
    fasta2peptides.py \\
        -i ${fasta} \\
        -o ${prefix} \\
        -minl ${min_length} \\
        -maxl ${max_length} \\
        -pepcol ${params.peptide_col_name} \\
        ${variant} \\
        ${wild_type} \\
        ${proteome}
    """

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def min_length = meta.mhc_class == "I" ? params.min_peptide_length_classI : params.min_peptide_length_classII
    def max_length = meta.mhc_class == "I" ? params.max_peptide_length_classI : params.max_peptide_length_classII
    """
    touch ${prefix}_length_${min_length}.tsv
    touch ${prefix}_length_${max_length}.tsv
    ${variants_tsv ? "touch ${prefix}.annotated.fasta" : ''}
    """
}
