process NETMHCIIPAN {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
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
    // Adjust for netMHCIIpan allele format (e.g. DRB1_0101, HLA-DPA10103-DPB10101, H-2-IAb)
    def alleles = meta.alleles_supported.tokenize(';')
                    .collect { allele ->
                        if (allele.contains('DRB')) {
                            // HLA-DRB1*01:01 -> DRB1_0101
                            allele.replace('*', '_').replace(':', '').replace('HLA-', '')
                        } else if (allele.startsWith('H2-') && allele.contains('/')) {
                            // mhcgnomes mouse class II canonical form is H2-<X>A*<Y>/<X>B*<Y>
                            // (locus X = A or E, haplotype Y = b/d/k/...). NetMHCIIpan wants H-2-I<X><Y>.
                            "H-2-I${allele[3]}${allele.substring(allele.lastIndexOf('*') + 1)}"
                        } else {
                            // HLA-DPA1*01:03/DPB1*04:01 -> HLA-DPA10103-DPB10401
                            allele.replace('*', '').replace(':', '').replace('/','-').replace('H2','H-2')
                        }
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
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_netmhciipan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhciipan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
