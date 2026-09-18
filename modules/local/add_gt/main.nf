process ADD_GT {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf)

    output:
    tuple val(meta), path("*.gt.vcf.gz"), path("*.gt.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args   = task.ext.args ?: '-t a -n c:0/1'
    def tumor  = meta.tumor_sample ?: ''
    """
    # pvacseq refuses VCFs without GT (Strelka emits none) and reads genotypes from the tumor
    # sample only, so setting every sample to 0/1 is enough here.
    tumor="${tumor}"
    if [ -n "\${tumor}" ] && ! bcftools query -l ${vcf} | grep -qx "\${tumor}"; then
        echo "ERROR: sample '\${tumor}' not found in ${vcf}." >&2
        exit 1
    fi
    if [ -z "\${tumor}" ] && [ "\$(bcftools query -l ${vcf} | wc -l)" -gt 1 ]; then
        echo "ERROR: ${vcf} has more than one sample; set tumor_sample in the samplesheet." >&2
        exit 1
    fi

    if bcftools view -h ${vcf} | grep -q '^##FORMAT=<ID=GT,'; then
        bcftools view ${vcf} -Oz -o ${prefix}.gt.vcf.gz
    else
        bcftools +setGT ${vcf} -Oz -o ${prefix}.gt.vcf.gz -- ${args}
    fi
    bcftools index -t ${prefix}.gt.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.gt.vcf.gz
    touch ${prefix}.gt.vcf.gz.tbi
    """
}
