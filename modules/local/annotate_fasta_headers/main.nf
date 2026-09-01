process ANNOTATE_FASTA_HEADERS {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    // The script is stdlib-only; reusing the sibling biopython image saves a container.
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/83/8372f6241b480332d91bc00a88ec8c72c8f7fcc9994177a5dd67a07007cd6e32/data' :
        'community.wave.seqera.io/library/biopython:1.85--6f761292fa9881b4' }"

    input:
    tuple val(meta), path(raw_fasta), path(vep_vcf)

    output:
    tuple val(meta), path("*.variant_peptides.annotated.fasta"), emit: fasta
    tuple val("${task.process}"), val('python'), eval("python3 --version | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    annotate_fasta_headers.py \\
        --vep-vcf ${vep_vcf} \\
        --in-fasta ${raw_fasta} \\
        --out-fasta ${prefix}.variant_peptides.annotated.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_peptides.annotated.fasta
    """
}
