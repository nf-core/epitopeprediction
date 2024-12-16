process NETMHCPAN {
    label 'process_single'
    tag "${meta.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file), path(software)

    output:
    tuple val(meta), path("*.xls"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II."
    }
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    // A*01:217 to HLA-A01:217 for meta.alleles: Add HLA- to the allele names and strip the *.
    def alleles    = meta.alleles.tokenize(';').collect { 'HLA-' + it.replace('*', '') }.join(',')

    """
    netmhcpan/netMHCpan \
        -p $peptide_file \
        -a $alleles \
        -xls \
        -xlsfile ${prefix}_predicted_netmhcpan.xls \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_netmhcpan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
