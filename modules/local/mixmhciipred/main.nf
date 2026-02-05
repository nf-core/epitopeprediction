process MIXMHCIIPRED {
    label 'process_single'
    tag "${meta.id}"

    // Container built on-the-fly via Wave from Dockerfile in this module directory.
    // NOTE: Do NOT add a container directive here - Wave will automatically use
    // the Dockerfile when wave.enabled = true.
    // MixMHCIIpred is NOT distributed due to license restrictions.
    // Requires: -profile wave (or wave.enabled = true in config)

    input:
    tuple val(meta), path(tsv)

    output:
    tuple val(meta), path("*_mixmhciipred.txt"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "II") {
        error "MIXMHCIIPRED only supports MHC class II. Use MIXMHCPRED for MHC class I."
    }
    def args   = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"
    // Convert alleles from mhcgnomes format to MixMHCIIpred format:
    // HLA-DRB1*03:01 -> DRB1_03_01
    // HLA-DPA1*01:03-DPB1*04:01 -> DPA1_01_03__DPB1_04_01
    def alleles = meta.alleles_supported.tokenize(';')
        .collect { allele ->
            allele.replace('HLA-','')
                  .replace('*','_')
                  .replace(':','_')
                  .replace('-','__')
        }
        .join(' ')
    """
    # MixMHCIIpred is pre-installed in container at /opt/MixMHC2pred
    # Use MixMHC2pred_unix for Linux systems
    /opt/MixMHC2pred/MixMHC2pred_unix \\
        -i $tsv \\
        -o ${prefix}_mixmhciipred.txt \\
        -a $alleles \\
        --no_context \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixmhciipred: "2.0.2"
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_mixmhciipred.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        mixmhciipred: "2.0.2"
    END_VERSIONS
    """
}
