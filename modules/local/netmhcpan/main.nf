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
    # TEMPORARY DIAGNOSTIC BLOCK -- to be removed before this PR is un-drafted.
    # All netMHCpan stdout is redirected to files: Nextflow only shows the TAIL of
    # .command.out, and the per-peptide table drowns the verdicts otherwise.
    # NB: no `| head` anywhere -- scripts run under `bash -ue -o pipefail`, so
    # SIGPIPE on the producer aborts the whole script with exit 141.
    sw=\$(readlink -f netmhcpan)
    in_abs=\$(readlink -f $tsv)
    sw_size=\$(du -sh --dereference netmhcpan | cut -f1)

    netmhcpan/netMHCpan -p $tsv -a $alleles -xls \
        -xlsfile ${prefix}_predicted_netmhcpan.xls $args > run1.log 2>&1 && r1=0 || r1=\$?

    if [ ! -s ${prefix}_predicted_netmhcpan.xls ]; then
        rm -rf /tmp/nmpA /tmp/nmpB
        mkdir -p /tmp/nmpA /tmp/nmpB

        ln -s "\$sw" /tmp/nmpA/nm
        ln -s "\$in_abs" /tmp/nmpA/in.tsv
        ( cd /tmp/nmpA && ./nm/netMHCpan -p in.tsv -a $alleles -xls -xlsfile probeA.out $args ) > probeA.log 2>&1 && rA=0 || rA=\$?

        cp -rL netmhcpan /tmp/nmpB/nm
        ln -s "\$in_abs" /tmp/nmpB/in.tsv
        ( cd /tmp/nmpB && ./nm/netMHCpan -p in.tsv -a $alleles -xls -xlsfile probeB.out $args ) > probeB.log 2>&1 && rB=0 || rB=\$?

        echo "===== VERDICTS ====="
        echo "PWD_LEN=\${#PWD}"
        echo "software_dir_size=\$sw_size"
        echo "symlink_target_len=\${#sw}"
        echo "run1(long cwd)      exit=\$r1 xls_bytes=\$(stat -c%s ${prefix}_predicted_netmhcpan.xls 2>/dev/null || echo 0)"
        echo "probeA(symlinked sw) exit=\$rA out_bytes=\$(stat -c%s /tmp/nmpA/probeA.out 2>/dev/null || echo 0)"
        echo "probeB(copied sw)    exit=\$rB out_bytes=\$(stat -c%s /tmp/nmpB/probeB.out 2>/dev/null || echo 0)"
        echo "probeA tail: \$(tail -2 probeA.log | tr '\\n' ' ')"
        echo "===================="
    fi

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
