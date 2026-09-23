process DOWNLOAD_REF_FASTA {
    tag "${meta.id}"
    label 'process_low'

    // conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://community-cr-prod.seqera.io/docker/registry/v2/blobs/sha256/a0/a01624095a85540784ea12ef530a030c33a34a7c51cccc16477dfba3a466d9a5/data'
        : 'community.wave.seqera.io/library/wget_gzip:174d767f72b71070'}"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path("${prefix}.fa"), emit: fasta
    tuple val("${task.process}"), val('wget'), eval("wget --version | head -n1 | sed 's/^GNU Wget //; s/ .*//'"), topic: versions, emit: versions_wget
    tuple val("${task.process}"), val('gzip'), eval("gzip --version | head -n1 | sed 's/^gzip //'"), topic: versions, emit: versions_gzip

    when:
    task.ext.when == null || task.ext.when

    script:
    // Fetched straight from Ensembl: vep_install's --AUTO f silently no-ops for many species.
    prefix = task.ext.prefix ?: "${species}.${assembly}"
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def sp_cap = species.capitalize()
    """
    # GRCh37 lives under a dedicated Ensembl FTP tree; everything else under the release root.
    if [ "${assembly}" = "GRCh37" ]; then
        base="https://ftp.ensembl.org/pub/grch37/release-${cache_version}"
    else
        base="https://ftp.ensembl.org/pub/release-${cache_version}"
    fi
    # Some genomes only ship a toplevel FASTA, so fall back to it.
    got=""
    for kind in primary_assembly toplevel; do
        url="\${base}/fasta/${species}/dna/${sp_cap}.${assembly}.dna.\${kind}.fa.gz"
        echo "Trying \${url}" >&2
        if wget ${args} -q -t 3 --timeout=60 -O ${prefix}.fa.gz "\${url}"; then got="\${kind}"; break; fi
        rm -f ${prefix}.fa.gz
    done
    if [ -z "\${got}" ]; then
        echo "ERROR: could not download a reference FASTA for ${species} ${assembly} (release ${cache_version}) from Ensembl." >&2
        exit 1
    fi
    echo "Downloaded \${got} FASTA" >&2

    gunzip ${args2} -f ${prefix}.fa.gz
    """

    stub:
    prefix = task.ext.prefix ?: "${species}.${assembly}"
    """
    touch ${prefix}.fa
    """
}
