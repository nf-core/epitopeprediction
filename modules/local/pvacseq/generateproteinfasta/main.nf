process PVACSEQ_GENERATEPROTEINFASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/pvactools:7.0.1--pyhdfd78af_0'
        : 'biocontainers/pvactools:7.0.1--pyhdfd78af_0'}"

    input:
    tuple val(meta), path(vcf), path(tbi), path(proximal_vcf), path(proximal_tbi)

    output:
    tuple val(meta), path("*.windows.fasta"), path("*.variants.tsv"), emit: fasta
    tuple val("${task.process}"), val('pvactools'), eval("pip show pvactools | grep '^Version:' | cut -d' ' -f2"), topic: versions, emit: versions_pvactools

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix = task.ext.prefix ?: "${meta.id}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def flank = params.mutation_flanking_aas
    def sample_arg = meta.tumor_sample ? "-s ${meta.tumor_sample}" : ''
    """
    # Each variant alone, then with its nearby variants folded in (assumed cis). Both sets are
    # kept, and stay in separate files so each mutant window is compared with the wild-type
    # window from the same run.
    pvacseq generate_protein_fasta ${vcf} ${flank} ${prefix}.1.windows.fasta ${sample_arg} ${args}
    pvacseq generate_protein_fasta ${vcf} ${flank} ${prefix}.2.windows.fasta ${sample_arg} ${args2} -p ${proximal_vcf}

    # pvacseq deletes the table its own FASTA ids are built from, so write it out here.
    pvacseq_variants_tsv.py \\
        --vep-vcf ${vcf} \\
        --output ${prefix}.variants.tsv \\
        ${sample_arg.replace('-s ', '--sample-name ')}
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.1.windows.fasta ${prefix}.2.windows.fasta ${prefix}.variants.tsv
    """
}
