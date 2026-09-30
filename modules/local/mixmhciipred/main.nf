process MIXMHCIIPRED {
    label 'process_single'
    tag "${meta.id}"

    // No container directive: the MixMHC2pred license prohibits redistribution, so Wave builds the
    // image on the fly from this module's Dockerfile (requires `-with-wave`).

    input:
    tuple val(meta), val(alleles_input), path(tsv)

    output:
    tuple val(meta), path("*_predicted_mixmhciipred.txt"), emit: predicted
    tuple val("${task.process}"), val('mixmhc2pred'), eval("sed -n 's/^# Output from MixMHC2pred (v\\(.*\\))/\\1/p' *_predicted_mixmhciipred.txt"), topic: versions, emit: versions_mixmhc2pred

    script:
    if (meta.mhc_class != "II") {
        error "MIXMHCIIPRED only supports MHC class II. Use MIXMHCPRED for MHC class I."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MixMHC2pred_unix \\
        -i $tsv \\
        -o ${prefix}_predicted_mixmhciipred.txt \\
        -a $alleles_input \\
        --no_context \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    echo '# Output from MixMHC2pred (v2.0.2)' > ${prefix}_predicted_mixmhciipred.txt
    """
}
