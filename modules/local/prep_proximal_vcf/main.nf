process PREP_PROXIMAL_VCF {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.proximal.vcf.gz"), path("*.proximal.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor  = meta.tumor_sample ?: ''
    """
    # pvacseq folds proximal variants into a window only when their HP phasing tag matches the main
    # variant's. Real phasing (HP already present) is kept; otherwise the same HP on every record
    # of the tumor sample treats all of them as cis.
    tumor="${tumor}"
    [ -n "\${tumor}" ] || tumor=\$(bcftools query -l ${vcf} | head -n 1)
    if bcftools view -h ${vcf} | grep -q '^##FORMAT=<ID=HP,'; then
        bcftools view -s "\${tumor}" ${vcf} -Oz -o ${prefix}.proximal.vcf.gz
    else
        bcftools view -s "\${tumor}" ${vcf} -Oz -o tumor.vcf.gz
        bcftools query -f '%CHROM\\t%POS\\t%REF\\t%ALT\\t1-1,1-2\\n' tumor.vcf.gz | bgzip > hp.tsv.gz
        tabix -s 1 -b 2 -e 2 hp.tsv.gz
        echo '##FORMAT=<ID=HP,Number=.,Type=String,Description="Read-backed phasing haplotype identifiers">' > hp.hdr
        bcftools annotate -a hp.tsv.gz -h hp.hdr -c CHROM,POS,REF,ALT,FMT/HP tumor.vcf.gz -Oz -o ${prefix}.proximal.vcf.gz
    fi
    bcftools index -t ${prefix}.proximal.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.proximal.vcf.gz
    touch ${prefix}.proximal.vcf.gz.tbi
    """
}
