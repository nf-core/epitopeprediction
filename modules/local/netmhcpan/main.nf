process NETMHCPAN {
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
    if (meta.mhc_class != "I") {
        error "NETMHCPAN only supports MHC class I. Use NETMHCIIPAN for MHC class II."
    }
    def args    = task.ext.args ?: ''
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def alleles = meta.alleles_supported.tokenize(';').collect { allele -> allele.replace('*', '').replace('H2','H-2') }.join(',')

    """
    # netMHCpan-4.2 builds internal paths from the current working directory into
    # fixed-size buffers. From a long cwd -- which Nextflow work dirs always are --
    # it overruns them and glibc aborts it ("*** buffer overflow detected ***").
    # The tcsh wrapper still exits 0, so the only symptom is a missing .xls.
    # Running from a short scratch dir avoids it; symlinks keep the ~195 MB tool
    # directory from being copied per task.
    scratch=\$(mktemp -d /tmp/netmhcpan.XXXXXX)
    ln -s "\$(readlink -f netmhcpan)" "\$scratch/nm"
    ln -s "\$(readlink -f $tsv)" "\$scratch/input.tsv"

    (
        cd "\$scratch"
        ./nm/netMHCpan \
            -p input.tsv \
            -a $alleles \
            -xls \
            -xlsfile ${prefix}_predicted_netmhcpan.xls \
            $args
    )

    # netMHCpan exits 0 even when it crashes, so check the output explicitly.
    if [ ! -s "\$scratch/${prefix}_predicted_netmhcpan.xls" ]; then
        echo "ERROR: netMHCpan produced no output for ${prefix}" >&2
        rm -rf "\$scratch"
        exit 1
    fi

    mv "\$scratch/${prefix}_predicted_netmhcpan.xls" .
    rm -rf "\$scratch"

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predicted_netmhcpan.xls

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        \$(cat netmhcpan/data/version | sed -s 's/ version/:/g')
    END_VERSIONS
    """
}
