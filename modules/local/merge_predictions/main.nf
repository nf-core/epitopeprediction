process MERGE_PREDICTIONS {
    label 'process_single'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    // Multi-arch Seqera Containers (mhcgnomes 3.32.0 + pyarrow 24.0.0):
    //   docker linux/amd64:        community.wave.seqera.io/library/mhcgnomes_pyarrow:6607e69dcab832f1
    //   docker linux/arm64:        community.wave.seqera.io/library/mhcgnomes_pyarrow:c31e3c4fb8f3dd4c
    //   singularity linux/amd64:   oras://community.wave.seqera.io/library/mhcgnomes_pyarrow:5069f7e652aac0d9
    //   singularity linux/arm64:   oras://community.wave.seqera.io/library/mhcgnomes_pyarrow:55b64b50976ff688
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'oras://community.wave.seqera.io/library/mhcgnomes_pyarrow:5069f7e652aac0d9' :
        'community.wave.seqera.io/library/mhcgnomes_pyarrow:6607e69dcab832f1' }"

    input:
    tuple val(meta), path(prediction_files), path(source_file)

    output:
    tuple val(meta), path("*.parquet") , emit: merged
    path "versions.yml"            , emit: versions

    script:
    template "merge_predictions.py"

    stub:
    def prefix     = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}_predictions.parquet

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version 2>&1 | sed 's/Python //g')
        pandas: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pandas').version)")
        mhcgnomes: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('mhcgnomes').version)")
        pyarrow: \$(python -c "import pkg_resources; print(pkg_resources.get_distribution('pyarrow').version)")
    END_VERSIONS
    """
}
