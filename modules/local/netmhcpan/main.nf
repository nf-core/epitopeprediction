process NETMHCPAN {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/de/de9c5fbcc5583f3c096617ef2c8f84c5e69b479cc5a5944f10d0e1d226779662/data' :
        'community.wave.seqera.io/library/bash_gawk_perl_tcsh:a941b4e9bd4b8805' }"

    input:
    tuple val(meta), val(alleles_input), path(tsv), path(software)

    output:
    tuple val(meta), path("*.xls"), emit: predicted
    tuple val("${task.process}"), val('netMHCpan'), eval("sed 's/.*version //' netmhcpan/data/version"), topic: versions

    script:
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II."
    }
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    // netMHCpan copies its install dir (NMHOME) and TMPDIR into fixed-size buffers (~95 chars) and aborts on long
    // work dir paths, so it is run through a short /tmp symlink with TMPDIR pointed there. See #341.
    """
    nm=\$(mktemp -d /tmp/nm.XXXXXX)
    trap 'rm -rf "\$nm"' EXIT
    ln -s "\$PWD/netmhcpan" "\$nm/netmhcpan"
    export TMPDIR="\$nm"

    "\$nm/netmhcpan/netMHCpan" \
        -p $tsv \
        -a $alleles_input \
        -xls \
        -xlsfile ${prefix}_predicted_netmhcpan.xls \
        $args
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_netmhcpan.xls
    """
}
