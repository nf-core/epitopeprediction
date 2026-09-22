process PREP_GERMLINE_CONTEXT {
    tag "${meta.id}"
    label 'process_single'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/0b/0b4d52ca9a56d07be3f78a12af654e5116f5112908dba277e6796fd9dfb83fe5/data'
        : 'community.wave.seqera.io/library/bcftools_htslib:1.23.1--9f08ec665533d64a'}"

    input:
    tuple val(meta), path(vcf), path(germline_vcf)
    path chr_map

    output:
    tuple val(meta), path("*.context.vcf.gz"), path("*.context.vcf.gz.tbi"), emit: vcf
    tuple val("${task.process}"), val('bcftools'), eval("bcftools --version | head -n1 | sed 's/^bcftools //'"), topic: versions, emit: versions_bcftools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: '-f PASS,.'
    def tumor = meta.tumor_sample ?: ''
    // pvacseq reads proximal variants within (flank + 1) * 4 bases of each somatic site; keeping
    // only those spares the second VEP run a whole germline call set.
    def window = (params.mutation_flanking_aas.toInteger() + 1) * 4
    """
    tumor="${tumor}"
    [ -n "\${tumor}" ] || tumor=\$(bcftools query -l ${vcf} | head -n 1)

    if [ "\$(bcftools query -l ${germline_vcf} | wc -l)" -ne 1 ]; then
        echo "ERROR: ${germline_vcf} must hold exactly one sample (the patient's normal)." >&2
        exit 1
    fi

    bcftools query -f '%CHROM\\t%POS\\n' ${vcf} \\
        | awk -v w=${window} -v OFS='\\t' '{ start = \$2 - w; if (start < 1) start = 1; print \$1, start, \$2 + w }' \\
        > windows.bed

    # both files must carry the same single sample to be concatenated
    bcftools view -s "\${tumor}" ${vcf} -Oz -o somatic.vcf.gz
    bcftools index -t somatic.vcf.gz

    # the somatic calls are already Ensembl-named at this point, so the germline ones must be too
    printf '%s\\n' "\${tumor}" > rename.txt
    bcftools view ${args} ${germline_vcf} -Ou \\
        | bcftools annotate --rename-chrs ${chr_map} -Ou \\
        | bcftools norm -m- -Oz -o germline.vcf.gz
    bcftools reheader -s rename.txt germline.vcf.gz -o germline.renamed.vcf.gz
    bcftools index -t germline.renamed.vcf.gz

    comm -12 <(bcftools query -f '%CHROM\\n' ${vcf} | sort -u) \\
             <(bcftools query -f '%CHROM\\n' germline.renamed.vcf.gz | sort -u) > shared_contigs.txt
    if [ ! -s shared_contigs.txt ]; then
        echo "ERROR: ${germline_vcf} shares no contig names with ${vcf}; is it the same genome build?" >&2
        exit 1
    fi

    bcftools view -R windows.bed germline.renamed.vcf.gz -Oz -o germline.near.vcf.gz
    bcftools index -t germline.near.vcf.gz
    kept=\$(bcftools view -H germline.near.vcf.gz | wc -l)
    echo "Germline records kept as context: \${kept}" >&2
    if [ "\${kept}" -eq 0 ]; then
        echo "WARNING: no germline record falls within a somatic window; the run proceeds without context." >&2
    fi

    bcftools concat -a somatic.vcf.gz germline.near.vcf.gz -Ou \\
        | bcftools sort -Oz -o ${prefix}.context.vcf.gz
    bcftools index -t ${prefix}.context.vcf.gz
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo | gzip > ${prefix}.context.vcf.gz
    touch ${prefix}.context.vcf.gz.tbi
    """
}
