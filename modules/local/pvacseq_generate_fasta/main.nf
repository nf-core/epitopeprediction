process PVACSEQ_GENERATE_FASTA {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pvactools:7.0.1--pyhdfd78af_0'
        : 'biocontainers/pvactools:7.0.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(proximal_vcf), path(proximal_tbi), val(min_length), val(max_length)

    output:
    tuple val(meta), path("*.len*.fasta"), path("*.variants.tsv"), emit: fasta
    tuple val("${task.process}"), val('pvactools'), eval("pip show pvactools | grep '^Version:' | cut -d' ' -f2"), topic: versions, emit: versions_pvactools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix      = task.ext.prefix ?: "${meta.id}"
    def args        = task.ext.args ?: ''
    def sample_name = meta.tumor_sample ?: ''
    def sample_arg  = sample_name ? "-s ${sample_name}" : ''
    def tsv_arg     = sample_name ? "--sample-name ${sample_name}" : ''
    """
    # As in `pvacseq run`: windows cut with k-1 flanking residues for peptide length k, so every
    # k-mer covers the mutation. Run 1 takes each variant alone, run 2 folds in nearby variants.
    for k in \$(seq ${min_length} ${max_length}); do
        pvacseq generate_protein_fasta ${vcf} \$((k - 1)) ${prefix}.len\${k}.1.fasta ${sample_arg} ${args}
        pvacseq generate_protein_fasta ${vcf} \$((k - 1)) ${prefix}.len\${k}.2.fasta ${sample_arg} ${args} -p ${proximal_vcf}
    done

    # pvacseq deletes the table its own FASTA ids are built from, so write it out here.
    pvacseq_variants_tsv.py \\
        --vep-vcf ${vcf} \\
        --output ${prefix}.variants.tsv \\
        ${tsv_arg}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.len${min_length}.1.fasta ${prefix}.len${max_length}.2.fasta ${prefix}.variants.tsv
    """
}
