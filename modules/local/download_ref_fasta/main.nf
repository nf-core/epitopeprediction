process DOWNLOAD_REF_FASTA {
    tag "${meta.id}"
    label 'process_low'

    conda "${moduleDir}/environment.yml"
    container "${workflow.containerEngine in ['singularity', 'apptainer'] && !task.ext.singularity_pull_docker_container
        ? 'https://depot.galaxyproject.org/singularity/samtools:1.21--h50ea8bc_0'
        : 'biocontainers/samtools:1.21--h50ea8bc_0'}"

    input:
    tuple val(meta), val(assembly), val(species), val(cache_version)

    output:
    tuple val(meta), path("${prefix}.fa")    , emit: fasta
    tuple val(meta), path("${prefix}.fa.fai"), emit: fai
    path "versions.yml"                      , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    // Fetch the reference FASTA straight from Ensembl rather than via vep_install's flaky
    // --AUTO f step, which silently no-ops for many species. Matches the cache build/version.
    prefix = task.ext.prefix ?: "${species}.${assembly}"
    def args = task.ext.args ?: ''
    """
    # GRCh37 lives under a dedicated Ensembl FTP tree; everything else under the release root.
    if [ "${assembly}" = "GRCh37" ]; then
        base="https://ftp.ensembl.org/pub/grch37/release-${cache_version}"
    else
        base="https://ftp.ensembl.org/pub/release-${cache_version}"
    fi
    sp="${species}"
    sp_cap="\${sp^}"   # homo_sapiens -> Homo_sapiens

    # Prefer the primary assembly (chromosomes + scaffolds, no patches/haplotypes); some
    # genomes only ship a toplevel FASTA, so fall back to it.
    got=""
    for kind in primary_assembly toplevel; do
        url="\${base}/fasta/${species}/dna/\${sp_cap}.${assembly}.dna.\${kind}.fa.gz"
        echo "Trying \${url}" >&2
        if wget ${args} -q -t 3 --timeout=60 -O ${prefix}.fa.gz "\${url}"; then got="\${kind}"; break; fi
        rm -f ${prefix}.fa.gz
    done
    if [ -z "\${got}" ]; then
        echo "ERROR: could not download a reference FASTA for ${species} ${assembly} (release ${cache_version}) from Ensembl." >&2
        exit 1
    fi
    echo "Downloaded \${got} FASTA" >&2

    # One-off decompression to a plain, indexed FASTA (the --ref_fasta contract shared with
    # the user-provided path, so PREP_VCF/VEP are identical across both).
    gunzip -f ${prefix}.fa.gz
    samtools faidx ${prefix}.fa

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """

    stub:
    prefix = task.ext.prefix ?: "${species}.${assembly}"
    """
    touch ${prefix}.fa ${prefix}.fa.fai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        samtools: \$(samtools --version | head -n1 | sed 's/^samtools //')
    END_VERSIONS
    """
}
