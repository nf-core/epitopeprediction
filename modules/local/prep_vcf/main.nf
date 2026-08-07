process PREP_VCF {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.prep.vcf.gz"), path("*.prep.vcf.gz.tbi"), emit: vcf
    path "versions.yml"                                              , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    # Build a contig rename map from the VCF's own ##contig headers (chr1->1, chrM->MT)
    # so records match the Ensembl-named VEP cache. Already-Ensembl VCFs map to themselves (no-op).
    bcftools view -h ${vcf} \\
        | awk -F'[<,=>]' '/^##contig/{for(i=1;i<=NF;i++) if(\$i=="ID"){c=\$(i+1); e=c; sub(/^chr/,"",e); if(e=="M")e="MT"; print c"\\t"e}}' \\
        > chr_map.txt

    # PASS-filter, rename contigs, split multiallelics and left-align/normalize against the reference.
    bcftools view -f PASS ${vcf} -Ou \\
        | bcftools annotate --rename-chrs chr_map.txt -Ou \\
        | bcftools norm -m- -f ${fasta} -Oz -o ${prefix}.prep.vcf.gz
    bcftools index -t ${prefix}.prep.vcf.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.prep.vcf.gz
    touch ${prefix}.prep.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^bcftools //')
    END_VERSIONS
    """
}
