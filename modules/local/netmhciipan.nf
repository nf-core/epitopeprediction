process NETMHCIIPAN {
    label 'process_single'
    tag "${meta.sample}"

    container 'ghcr.io/jonasscheid/epitopeprediction-2:netmhc'

    input:
    tuple val(meta), path(peptide_file), path(software)

    output:
    tuple val(meta), path("*.xls"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCIIPAN for MHC class II."
    }
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: "${meta.id}"
    // Adjust for netMHCIIpan allele format (e.g. DRB1_0101, HLA-DPA10103-DPB10101)
    def alleles = meta.alleles_supported.tokenize(';')
                    .collect {
                        it.contains('DRB') ?
                            it.replace('*', '_').replace(':', '').replace('HLA-', '') :
                            it.replace('*', '').replace(':', '').replace('/','-')
                    }.join(',')

    """
    netmhciipan/netMHCIIpan \
        -f $peptide_file \
		-inptype 1 \
        -a $alleles \
        -xls \
        -xlsfile ${prefix}_predicted_netmhciipan.xls \
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhciipan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """

    stub:
    def args       = task.ext.args ?: ''
    def prefix     = task.ext.prefix ?: meta.sample
    """
    touch ${prefix}_predicted_netmhciipan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhciipan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
