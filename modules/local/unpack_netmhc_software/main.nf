/*
* Copy non-free software provided by the user into the working directory
*/
process UNPACK_NETMHC_SOFTWARE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container
        ? 'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img'
        : 'docker.io/biocontainers/biocontainers:v1.2.0_cv2'}"

    input:
    tuple val(toolname), val(toolversion), path(tooltarball), val(toolbinaryname)

    output:
    path "${toolname}", emit: nonfree_tools
    tuple val("${task.process}"), val(toolname), eval("sed 's/.*version //' ${toolname}/data/version"), topic: versions, emit: versions_netmhc

    when:
    task.ext.when == null || task.ext.when

    script:
    def expected_dir = "${toolbinaryname}-${toolversion}"
    """
    fail() {
        echo "Invalid ${toolname} tarball '${tooltarball}': \$1" >&2
        echo "Please provide an original ${toolbinaryname} ${toolversion} Linux tarball, any sub-release (e.g. ${toolbinaryname}-${toolversion}b.Linux.tar.gz) is accepted." >&2
        exit 1
    }

    [ -f "${tooltarball}" ] || fail "not a regular file."
    tar -tzf "${tooltarball}" > contents.txt 2> /dev/null || fail "not a readable .tar.gz archive (incomplete download or a saved web page?)."

    top_dirs="\$(awk -F/ '\$1 !~ /^\\./ { print \$1 }' contents.txt | sort -u | tr '\\n' ' ')"
    [ "\$top_dirs" = "${expected_dir} " ] || fail "expected the top-level folder '${expected_dir}', found: \${top_dirs% }"
    grep -q "^${expected_dir}/Linux_x86_64/bin/" contents.txt || fail "no Linux_x86_64 binaries found, only the Linux tarball is supported."
    grep -qx "${expected_dir}/data/version" contents.txt || fail "no data/version file found."

    mkdir -v "${toolname}"
    tar -C "${toolname}" --strip-components 1 -x -f "${tooltarball}"

    version="\$(sed -n 's/^.* version //p' "${toolname}/data/version")"
    case "\$version" in
        ${toolversion} | ${toolversion}[a-z]*) echo "Found ${toolbinaryname} version \$version" ;;
        *) fail "data/version reports '\$version', expected ${toolversion} or one of its sub-releases." ;;
    esac

    #
    # MODIFY THE NETMHC WRAPPER SCRIPT ACCORDING TO INSTALL INSTRUCTIONS
    # Substitution 1: We install tcsh via conda, thus /bin/tcsh won't work
    # Substitution 2: We want temp files to be written to /tmp if TMPDIR is not set
    # Substitution 3: NMHOME should be the folder in which the tcsh script itself resides
    #
    sed -i.bak \
        -e 's_bin/tcsh.*\$_usr/bin/env tcsh_' \
        -e "s_/scratch_/tmp_" \
        -e "s_setenv[[:space:]]NMHOME.*_setenv NMHOME \\`realpath -s \\\$0 | sed -r 's/[^/]+\$//'\\`_ " "${toolname}/${toolbinaryname}"
    """

    stub:
    """
    mkdir -p "${toolname}/data"
    echo "${toolbinaryname} version ${toolversion}" > "${toolname}/data/version"
    """
}
