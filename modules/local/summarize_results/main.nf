process SUMMARIZE_RESULTS {
    label 'process_low'
    tag "${meta.id}"

    // conda "${moduleDir}/environment.yml"
    // Multi-arch Seqera Containers (mhcgnomes 3.32.0 + pyarrow 24.0.0):
    //   docker linux/amd64:        community.wave.seqera.io/library/mhcgnomes_pyarrow:6607e69dcab832f1
    //   docker linux/arm64:        community.wave.seqera.io/library/mhcgnomes_pyarrow:c31e3c4fb8f3dd4c
    //   singularity linux/amd64:   https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7df37ea371a12296bc6064ad80c0bd27152061b8ca3f0a6ad3dd924deb813f0/data
    //   singularity linux/arm64:   https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/2d/2d0f348e08592099a98cd5d85ae7a2024e023dd4a95b593b356781ab4d7df6f5/data
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/e7/e7df37ea371a12296bc6064ad80c0bd27152061b8ca3f0a6ad3dd924deb813f0/data' :
        'community.wave.seqera.io/library/mhcgnomes_pyarrow:6607e69dcab832f1' }"

    input:
    tuple val(meta), path(parquet)

    output:
    tuple val(meta), path("*.tsv")      , emit: tsv
    tuple val(meta), path("*_mqc.json") , emit: json
    path "versions.yml"                 , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args    = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    summarize_results.py \\
        --input . \\
        --prefix ${prefix} \\
        --peptide_col_name ${params.peptide_col_name} \\
        $args


    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"

    """
    touch ${prefix}.tsv
    touch ${prefix}_mqc.json

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        python: \$(python --version | sed 's/Python //g')
        pandas: \$(python -c "import pandas; print(pandas.__version__)")
        pyarrow: \$(python -c "import pyarrow; print(pyarrow.__version__)")
    END_VERSIONS
    """
}
