process ANNOTATE_FASTA_HEADERS {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    // annotate_fasta_headers.py is stdlib-only; reuse the sibling biopython image so the
    // variant path needs one fewer container. This is the single VCF-join site: it enriches
    // the pvacseq WT/MT windows with VEP provenance so VARIANT_FASTA2PEPTIDES needs no VCF.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(raw_fasta), path(vep_vcf)

    output:
    tuple val(meta), path("*.variant_peptides.annotated.fasta"), emit: fasta
    path "versions.yml"                                        , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_fasta_headers.py \\
        --vep-vcf ${vep_vcf} \\
        --in-fasta ${raw_fasta} \\
        --out-fasta ${prefix}.variant_peptides.annotated.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_peptides.annotated.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python3 --version 2>&1 | cut -d' ' -f2)
    END_VERSIONS
    """
}
