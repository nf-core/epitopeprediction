process PVACSEQ_GENERATE_FASTA {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pvactools:7.0.1--pyhdfd78af_0'
        : 'biocontainers/pvactools:7.0.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(proximal_vcf), path(proximal_tbi)

    output:
    tuple val(meta), path("*.raw.fasta"), path(vcf), emit: fasta
    tuple val("${task.process}"), val('pvactools'), eval("pip show pvactools | grep '^Version:' | cut -d' ' -f2"), topic: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    def args       = task.ext.args ?: ''
    def flank      = params.mutation_flanking_aas
    // -s picks the tumor column on multi-sample VCFs; single-sample VCFs leave it unset.
    def sample_arg = meta.tumor_sample ? "-s ${meta.tumor_sample}" : ''
    """
    # Each variant alone, then with its proximal variants folded in (assumed cis); keep both so
    # single- and multi-variant peptides are generated, dropping windows that came out identical.
    pvacseq generate_protein_fasta ${vcf} ${flank} single.fasta ${sample_arg} ${args}
    pvacseq generate_protein_fasta ${vcf} ${flank} proximal.fasta ${sample_arg} ${args} -p ${proximal_vcf}
    awk 'function emit() { if (header != "" && !seen[header SUBSEP seq]++) print header seq }
         /^>/ { emit(); header = \$0; seq = ""; next }
         { seq = seq "\\n" \$0 }
         END { emit() }' single.fasta proximal.fasta > ${prefix}.raw.fasta
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.raw.fasta
    """
}
