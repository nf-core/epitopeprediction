process PREP_VCF {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/3c/3cc9ab3025e57f15ff3d58e45dd31c746d0801b8f2456fbf39eb7db8368de4e8/data'
        : 'community.wave.seqera.io/library/bcftools_htslib_vatools:ef4839f79b3f57f7'}"

    input:
    tuple val(meta), path(vcf)
    tuple val(meta2), path(fasta)
    tuple val(meta3), path(fai)

    output:
    tuple val(meta), path("*.prep.vcf.gz"), path("*.prep.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions
    tuple val("${task.process}"), val('vatools'), eval("pip show vatools | sed -n 's/^Version: //p'"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def tumor  = meta.tumor_sample ?: ''
    """
    # pvacseq refuses VCFs without GT (Strelka emits none). Add GT=0/1 for the tumor sample with vatools,
    # as pVACtools recommends (https://pvactools.readthedocs.io/en/latest/pvacseq/input_file_prep/gt.html).
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
        elif ! bcftools query -l ${vcf} | grep -qx "\${tumor}"; then
            echo "ERROR: sample '\${tumor}' not found in ${vcf}." >&2
            exit 1
        fi
        vcf-genotype-annotator ${vcf} "\${tumor}" 0/1 -o ${prefix}.gt.vcf
        input_vcf=${prefix}.gt.vcf
    fi

    # Rename contigs to Ensembl style (chr1->1, chrM->MT) so records match the VEP cache; names come
    # from the ##contig headers, or from the records when a VCF has none.
    names=\$(bcftools view -h ${vcf} | awk -F'[<,=>]' '/^##contig/ { for (i = 1; i <= NF; i++) if (\$i == "ID") print \$(i + 1) }')
    [ -n "\${names}" ] || names=\$(bcftools query -f '%CHROM\\n' ${vcf} | sort -u)
    echo "\${names}" | awk '{ ensembl = \$1; sub(/^chr/, "", ensembl); if (ensembl == "M") ensembl = "MT"; print \$1 "\\t" ensembl }' > chr_map.txt

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
