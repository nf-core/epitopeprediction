process NETMHCPAN {
    label 'process_single'
    tag "${meta.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file), path(nonfree_tools)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II, or adjust the samplesheet accordingly."
    }
    def prefix = peptide_file.baseName

    """
    """
    stub:
    """
    touch ${meta.sample}_predicted_netmhcpan.tsv
    touch versions.yml
    """
}
