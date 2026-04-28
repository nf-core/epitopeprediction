process MIXMHCPRED {
    label 'process_single'
    tag "${meta.id}"

    // Container built on-the-fly via Wave from Dockerfile in this module directory.
    // NOTE: Do NOT add a container directive here - Wave will automatically use
    // the Dockerfile when wave.enabled = true.
    // MixMHCpred is NOT distributed due to license restrictions.
    // Requires: -profile wave (or wave.enabled = true in config)

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*_mixmhcpred.txt"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "I") {
        error "MIXMHCPRED only supports MHC class I."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Convert alleles from mhcgnomes format (HLA-A*01:01) to MixMHCpred format (A0101)
    def alleles = meta.alleles_supported.tokenize(';')
        .collect { it.replace('HLA-','').replace('*','').replace(':','') }
        .join(',')
    """
    # MixMHCpred is pre-installed in container at /opt/MixMHCpred
    /opt/MixMHCpred/MixMHCpred \\
        -i $tsv \\
        -o ${prefix}_mixmhcpred.txt \\
        -a $alleles \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixmhcpred: \$(/opt/MixMHCpred/MixMHCpred -h 2>&1 | head -1 | grep -oE '[0-9]+\\.[0-9.]+')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mixmhcpred.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixmhcpred: "3.0"
    END_VERSIONS
    """
}
