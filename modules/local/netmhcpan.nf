process NETMHCPAN {
    label 'process_single'
    tag "${meta.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file), path(software)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II."
    }
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample

    """
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_netmhcpan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhcpan \$(cat data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
