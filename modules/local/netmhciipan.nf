process NETMHCIIPAN {
    label 'process_single'
    tag "${metadata.sample}"


    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(metadata), path(peptide_file)

    output:
    tuple val(metadata), path("*.tsv"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (metadata.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCPAN for MHC class I, or adjust the samplesheet accordingly."
    }

    """
    """

    stub:
    """
    touch ${metadata.sample}_predicted_netmhciipan.tsv
    touch versions.yml
    """
}
