process MIXMHCPRED {
    label 'process_single'
    tag "${meta.id}"

    // No container directive: the MixMHCpred license prohibits redistribution, so Wave builds the
    // image on the fly from this module's Dockerfile (requires `-with-wave`).

    input:
    tuple val(meta), val(alleles_input), path(tsv)

    output:
    tuple val(meta), path("*_predicted_mixmhcpred.txt"), emit: predicted
    tuple val("${task.process}"), val('mixmhcpred'), eval("MixMHCpred -h | head -1 | sed 's/MixMHCpred//'"), topic: versions, emit: versions_mixmhcpred

    script:
    if (meta.mhc_class != "I") {
        error "MIXMHCPRED only supports MHC class I. Use MIXMHCIIPRED for MHC class II."
    }
    if (!task.container && !workflow.wave?.enabled) {
        log.warn1("MIXMHCPRED has no public container: add `-with-wave` to build it on the fly, or set your own container for the process.")
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    MixMHCpred \\
        -i $tsv \\
        -o ${prefix}_predicted_mixmhcpred.txt \\
        -a $alleles_input \\
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_mixmhcpred.txt
    """
}
