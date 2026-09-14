process PREP_VCF {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.prep.vcf.gz"), path("*.prep.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor  = meta.tumor_sample ?: ''
    """
    # pvacseq refuses VCFs without GT (Strelka emits none): tumor sample gets 0/1, the rest stay missing.
    if bcftools view -h ${vcf} | grep -q '^##FORMAT=<ID=GT,'; then
        input_vcf=${vcf}
    else
        tumor="${tumor}"
        if [ -z "\${tumor}" ]; then
            if [ "\$(bcftools query -l ${vcf} | wc -l)" -ne 1 ]; then
                echo "ERROR: ${vcf} has no GT field and more than one sample; set tumor_sample in the samplesheet." >&2
                exit 1
            fi
            tumor=\$(bcftools query -l ${vcf})
        fi
        mkdir gt
        bcftools view ${vcf} -Oz -o gt/input.vcf.gz && bcftools index -t gt/input.vcf.gz
        bcftools view -s "\${tumor}" gt/input.vcf.gz -Ou \\
            | bcftools +setGT -Oz -o gt/tumor.vcf.gz -- -t a -n c:0/1
        bcftools index -t gt/tumor.vcf.gz
        bcftools annotate -a gt/tumor.vcf.gz -c FMT/GT gt/input.vcf.gz -Oz -o gt/with_gt.vcf.gz
        input_vcf=gt/with_gt.vcf.gz
    fi

    # Rename map from the VCF's own ##contig headers (chr1->1, chrM->MT) so records match
    # the Ensembl-named VEP cache. Already-Ensembl VCFs map to themselves.
    bcftools view -h ${vcf} | awk -F'[<,=>]' '
        /^##contig/ {
            for (i = 1; i <= NF; i++) if (\$i == "ID") name = \$(i + 1)
            ensembl = name
            sub(/^chr/, "", ensembl)
            if (ensembl == "M") ensembl = "MT"
            print name "\\t" ensembl
        }' > chr_map.txt

    bcftools view -f PASS \${input_vcf} -Ou \\
        | bcftools annotate --rename-chrs chr_map.txt -Ou \\
        | bcftools norm -m- -f ${fasta} -Oz -o ${prefix}.prep.vcf.gz
    bcftools index -t ${prefix}.prep.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.prep.vcf.gz
    touch ${prefix}.prep.vcf.gz.tbi
    """
}
