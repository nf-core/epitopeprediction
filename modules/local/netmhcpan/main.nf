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
    # netMHCpan-4.2 aborts with "*** buffer overflow detected ***" after printing
    # "# Predict", exits 0 anyway (tcsh wrapper), so Nextflow only reports the
    # missing *.xls. The probes below narrow down what triggers the overflow.
    # NB: no `| head` anywhere below -- Nextflow runs this with `bash -ue -o pipefail`,
    # so SIGPIPE on the producer aborts the whole script with exit 141.
    echo "### env"
    echo "PWD_LEN=\${#PWD}"
    echo "TMPDIR=\${TMPDIR:-unset}"
    ldd --version > ldd.txt 2>&1 || true
    sed -n '1p' ldd.txt
    echo "stack_limit=\$(ulimit -s)"
    echo "### input"
    echo "alleles=$alleles"
    echo "n_peptides=\$(wc -l < $tsv)"
    echo "max_peptide_len=\$(awk '{ if (length(\$0) > m) m = length(\$0) } END { print m }' $tsv)"
    sed -n '1,3p' $tsv | cat -A
    tail -2 $tsv | cat -A

    echo "### run 1: as-is"
    netmhcpan/netMHCpan \
        -p $tsv \
        -a $alleles \
        -xls \
        -xlsfile ${prefix}_predicted_netmhcpan.xls \
        $args && rc=0 || rc=\$?
    echo "run1_exit=\$rc"
    ls -l ${prefix}_predicted_netmhcpan.xls 2>&1 || echo "run1: NO XLS"

    if [ ! -s ${prefix}_predicted_netmhcpan.xls ]; then
        echo "software_dir_size=\$(du -sh --dereference netmhcpan | cut -f1)"

        echo "### probe A: short cwd, software + input as SYMLINKS (cheap)"
        rm -rf /tmp/nmpA && mkdir -p /tmp/nmpA
        ln -s "\$(readlink -f netmhcpan)" /tmp/nmpA/nm
        ln -s "\$(readlink -f $tsv)" /tmp/nmpA/in.tsv
        ( cd /tmp/nmpA \
            && ./nm/netMHCpan -p in.tsv -a $alleles -xls -xlsfile probeA.out $args && rc=0 || rc=\$? \
            ; echo "probeA_exit=\$rc" \
            ; ls -l /tmp/nmpA/probeA.out 2>&1 || echo "probeA: NO OUTPUT" )

        echo "### probe B: short cwd, software COPIED, input symlinked"
        rm -rf /tmp/nmpB && mkdir -p /tmp/nmpB
        cp -rL netmhcpan /tmp/nmpB/nm
        ln -s "\$(readlink -f $tsv)" /tmp/nmpB/in.tsv
        ( cd /tmp/nmpB \
            && ./nm/netMHCpan -p in.tsv -a $alleles -xls -xlsfile probeB.out $args && rc=0 || rc=\$? \
            ; echo "probeB_exit=\$rc" \
            ; ls -l /tmp/nmpB/probeB.out 2>&1 || echo "probeB: NO OUTPUT" )
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
