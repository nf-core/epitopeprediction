process NETMHCIIPAN {
    label 'process_single'
    tag "${meta.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file), path(software)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCPAN for MHC class I."
    }
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample

    """
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_netmhciipan.tsv

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        netmhciipan \$(cat data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
