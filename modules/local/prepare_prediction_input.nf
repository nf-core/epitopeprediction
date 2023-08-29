//and also do the allele name

process PREPARE_PREDICITION_INPUT {
    label 'process_single'
    tag "${metadata.sample}"

    input:
    tuple val(metadata), path(prediction_files), path(metadata_file)

    output:
    tuple val(metadata), path("*.tsv"), emit: merged
    path "versions.yml", emit: versions

    script:
    """
    """

    stub:
    """
    """
}
