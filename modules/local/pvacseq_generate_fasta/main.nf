process PVACSEQ_GENERATE_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pvactools:7.0.1--pyhdfd78af_0'
        : 'biocontainers/pvactools:7.0.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi)

    output:
    tuple val(meta), path("*.variant_peptides.raw.fasta"), path(vcf), emit: fasta
    path "versions.yml"                                             , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args ?: ''
    def flank      = params.flank
    def downstream = params.downstream
    // Select the tumor sample's genotypes on multi-sample (matched tumor/normal) VCFs.
    // Optional: single-sample tumor-only VCFs leave meta.tumor_sample unset -> no -s.
    def sample_arg = meta.tumor_sample ? "-s ${meta.tumor_sample}" : ''
    """
    pvacseq generate_protein_fasta \\
        ${vcf} \\
        ${flank} \\
        ${prefix}.variant_peptides.raw.fasta \\
        -d ${downstream} \\
        --pass-only \\
        ${sample_arg} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pip show pvactools 2>/dev/null | awk '/^Version:/{print \$2}')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.variant_peptides.raw.fasta

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        pvactools: \$(pip show pvactools 2>/dev/null | awk '/^Version:/{print \$2}')
    END_VERSIONS
    """
}
