process NETMHCIIPAN {
    label 'process_single'
    tag "${meta.sample}"


    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file)

    output:
    tuple val(meta), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCPAN for MHC class I, or adjust the samplesheet accordingly."
    }

    """
    """

    stub:
    """
    touch ${meta.sample}_predicted_netmhciipan.tsv
    touch versions.yml
    """
}
