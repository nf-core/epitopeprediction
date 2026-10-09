/*
* Copy non-free software provided by the user into the working directory
*/
process UNPACK_NETMHC_SOFTWARE {
    label 'process_single'

    conda "${moduleDir}/environment.yml"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://containers.biocontainers.pro/s3/SingImgsRepo/biocontainers/v1.2.0_cv1/biocontainers_v1.2.0_cv1.img' :
        'docker.io/biocontainers/biocontainers:v1.2.0_cv2' }"

    input:
    tuple val(toolname), val(toolversion), path(tooltarball), val(toolbinaryname)

    output:
    path "${toolname}", emit: nonfree_tools
    tuple val("${task.process}"), val(toolname), val(toolversion), topic: versions, emit: versions_netmhc

    when:
    task.ext.when == null || task.ext.when

    script:
    """
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
    """

    stub:
    """
    mkdir "${toolname}"
    """
}
