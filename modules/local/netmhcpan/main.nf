process NETMHCPAN {
    label 'process_single'
    tag "${meta.sample}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/bash_gawk_perl_tcsh:59b949a8a9f0816b' :
        'community.wave.seqera.io/library/bash_gawk_perl_tcsh:a941b4e9bd4b8805' }"

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
    def alleles    = meta.alleles_supported.tokenize(';').collect { it.replace('*', '') }.join(',')

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
