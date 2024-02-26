/*
* Copy non-free software provided by the user into the working directory
*/
process EXTERNAL_TOOLS_IMPORT {
    label 'process_low'

    conda "conda-forge::coreutils=9.1"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv2' }"

    input:
    tuple val(toolname), val(toolversion), val(toolchecksum), path(tooltarball), file(datatarball), val(datachecksum), val(toolbinaryname)

    output:
    path "${toolname}", emit: nonfree_tools
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    """
    #
    # CHECK IF THE PROVIDED SOFTWARE TARBALL IS A REGULAR FILES
    #
    if [ ! -f "$tooltarball" ]; then
        echo "Path specified for ${toolname} does not point to a regular file. Please specify a path to the original tool tarball." >&2
        exit 1
    fi

    #
    # VALIDATE THE CHECKSUM OF THE PROVIDED SOFTWARE TARBALL
    #
    checksum="\$(md5sum "$tooltarball" | cut -f1 -d' ')"
    echo "\$checksum"
    if [ "\$checksum" != "${toolchecksum}" ]; then
        echo "Checksum error for $toolname. Please make sure to provide the original tarball for $toolname version $toolversion" >&2
        exit 2
    fi

    #
    # UNPACK THE PROVIDED SOFTWARE TARBALL
    #
    mkdir -v "${toolname}"
    tar -C "${toolname}" --strip-components 1 -x -f "$tooltarball"

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

    # MODIFY perl location in perl script for netmhcIIpan
    if [ "$toolname" == "netmhciipan" ]; then
        sed -i.bak \
        -e 's_bin/perl.*\$_local/bin/perl_' "${toolname}/NetMHCIIpan-${toolversion}.pl"
    fi
    #
    # VALIDATE THE CHECKSUM OF THE DOWNLOADED MODEL DATA
    #
    checksum="\$(md5sum "$datatarball" | cut -f1 -d' ')"
    if [ "\$checksum" != "${datachecksum}" ]; then
        echo "A checksum mismatch occurred when checking the data file for ${toolname}." >&2
        exit 3
    fi

    #
    # UNPACK THE DOWNLOADED MODEL DATA
    #
    tar -C "${toolname}" -v -x -f "$datatarball"

    #
    # CREATE VERSION FILE
    #
    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ${toolname}: ${toolversion}
    END_VERSIONS
    """
}
