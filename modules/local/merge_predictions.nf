process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${metadata.sample}"

    conda "bioconda::mhcgnomes=1.8.4"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/mhcgnomes:1.8.4--pyh7cba7a3_0' :
        'quay.io/biocontainers/mhcgnomes:1.8.4--pyh7cba7a3_0' }"

    input:
    tuple val(metadata), path(prediction_files), path(metadata_file)

    output:
    tuple val(metadata), path("*.tsv"), emit: merged
    path "versions.yml", emit: versions

    script:
    def output = prediction_files.first().baseName.split("_").dropRight(2).join("_")
    def min_length = (metadata.mhc_class == "I") ? params.min_peptide_length_mhc_I : params.min_peptide_length_mhc_II
    def max_length = (metadata.mhc_class == "I") ? params.max_peptide_length_mhc_I : params.max_peptide_length_mhc_II

    def syfpeithi_threshold = params.syfpeithi_threshold ? "--syfpeithi_threshold ${params.syfpeithi_threshold}" : ""
    def mhcflurry_threshold = params.mhcflurry_threshold ? "--mhcflurry_threshold ${params.mhcflurry_threshold}" : ""
    def mhcnuggets_threshold = params.mhcnuggets_threshold ? "--mhcnuggets_threshold ${params.mhcnuggets_threshold}" : ""
    def netmhcpan_threshold = params.netmhcpan_threshold ? "--netmhcpan_threshold ${params.netmhcpan_threshold}" : ""
    def netmhciipan_threshold = params.netmhciipan_threshold ? "--netmhciipan_threshold ${params.netmhciipan_threshold}" : ""

    """
    merge_binding_predictions.py \
        --input ${prediction_files} \
        --input_metadata ${metadata_file} \
        --output ${output}_merged.tsv \
        --min_peptide_length ${min_length} \
        --max_peptide_length ${max_length} \
        $syfpeithi_threshold \
        $mhcflurry_threshold \
        $mhcnuggets_threshold \
        $netmhcpan_threshold \
        $netmhciipan_threshold


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python \$(python --version | sed 's/Python //g')
        mhcgnomes \$(python -c "from mhcgnomes import version; print(version.__version__)"   )
    END_VERSIONS
    """

    stub:
    """
    touch merged_prediction.tsv
    touch versions.yml
    """
}
