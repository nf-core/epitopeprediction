process DOWNLOAD_VEP_CACHE {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path(prefix), emit: cache

    when:
    task.ext.when == null || task.ext.when

    script:
    // Fetched straight from Ensembl over HTTPS: vep_install lists caches over FTP and, when that
    // fails, silently falls back to a six-entry species list.
    prefix = task.ext.prefix ?: 'vep_cache'
    def args     = task.ext.args ?: ''
    def base     = "https://ftp.ensembl.org/pub/release-${cache_version}/variation/indexed_vep_cache"
    def filename = "${species}_vep_${cache_version}_${assembly}.tar.gz"
    """
    wget ${args} -q -T 60 -O CHECKSUMS "${base}/CHECKSUMS"
    if ! grep -q " ${filename}\$" CHECKSUMS; then
        echo "ERROR: ${filename} not found at ${base} (check --vep_species/--vep_genome/--vep_cache_version)" >&2
        exit 1
    fi

    wget ${args} -q -T 60 -O ${filename} "${base}/${filename}"

    # CHECKSUMS holds BSD 'sum' output: <checksum> <1K blocks> <file>
    expected=\$(grep " ${filename}\$" CHECKSUMS | awk '{print \$1, \$2}')
    actual=\$(sum ${filename} | awk '{print \$1, \$2}')
    if [ "\${expected}" != "\${actual}" ]; then
        echo "ERROR: checksum mismatch for ${filename}: expected '\${expected}', got '\${actual}'" >&2
        exit 1
    fi

    mkdir ${prefix}
    tar -xzf ${filename} -C ${prefix}
    rm ${filename}
    """

    stub:
    prefix = task.ext.prefix ?: 'vep_cache'
    """
    mkdir -p ${prefix}/${species}/${cache_version}_${assembly}
    """
}
