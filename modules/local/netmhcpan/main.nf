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
        echo "### run 2: single peptide, same directory"
        sed -n '1p' $tsv > one_peptide.txt
        netmhcpan/netMHCpan -p one_peptide.txt -a $alleles -xls -xlsfile probe_one.out $args && rc=0 || rc=\$?
        echo "run2_exit=\$rc"
        ls -l probe_one.out 2>&1 || echo "run2: NO OUTPUT"

        echo "### run 3: full input, short path (/tmp/nmp)"
        rm -rf /tmp/nmp && mkdir -p /tmp/nmp
        cp -rL netmhcpan /tmp/nmp/nm
        cp $tsv /tmp/nmp/in.tsv
        ( cd /tmp/nmp \
            && ./nm/netMHCpan -p in.tsv -a $alleles -xls -xlsfile probe_short.out $args && rc=0 || rc=\$? \
            ; echo "run3_exit=\$rc" \
            ; ls -l /tmp/nmp/probe_short.out 2>&1 || echo "run3: NO OUTPUT" )

        echo "### run 4: no -xls, stdout only"
        netmhcpan/netMHCpan -p $tsv -a $alleles $args > probe_stdout.txt 2>&1 && rc=0 || rc=\$?
        echo "run4_exit=\$rc"
        echo "run4_stdout_lines=\$(wc -l < probe_stdout.txt)"
        tail -5 probe_stdout.txt
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
