process NETMHCIIPAN {
    label 'process_single'
    tag "${meta.id}"

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de9c5fbcc5583f3c096617ef2c8f84c5e69b479cc5a5944f10d0e1d226779662/data' :
        'community.wave.seqera.io/library/bash_gawk_perl_tcsh:a941b4e9bd4b8805' }"

    input:
    tuple val(meta), path(tsv), path(software)

    output:
    tuple val(meta), path("*.xls"), emit: predicted
    path "versions.yml", emit: versions

    script:
    if (meta.mhc_class != "II") {
        error "NETMHCIIPAN only supports MHC class II. Use NETMHCPAN for MHC class I."
    }
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    // Adjust for netMHCIIpan allele format (e.g. DRB1_0101, HLA-DPA10103-DPB10101)
    def alleles = meta.alleles_supported.tokenize(';')
                    .collect {
                        it.contains('DRB') ?
                            it.replace('*', '_').replace(':', '').replace('HLA-', '') :
                            it.replace('*', '').replace(':', '').replace('/','-').replace('H2','H-2')
                    }.join(',')

    """
    netmhciipan/netMHCIIpan \
        -f $tsv \
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
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_netmhciipan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhciipan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
