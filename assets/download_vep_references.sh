#!/usr/bin/env bash
#
# Downloads the VEP offline cache (--vep_cache) and the reference genome FASTA + .fai
# (--ref_fasta) that the variant path needs, ~20 GB / ~26 GB extracted. Alternatively pass
# --vep_download_cache to fetch both in-pipeline. Needs curl, tar, gunzip and docker
# (docker only to samtools-faidx the FASTA).
#
# Usage:  ./download_vep_references.sh [REFDIR]        (default: <repo>/references)
#         SPECIES=mus_musculus ASSEMBLY=GRCm39 RELEASE=110 ./download_vep_references.sh
#         KEEP_ARCHIVES=1 ./download_vep_references.sh
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REFDIR="${1:-$(cd "$HERE/.." && pwd)/references}"

SPECIES="${SPECIES:-homo_sapiens}"
ASSEMBLY="${ASSEMBLY:-GRCh38}"
RELEASE="${RELEASE:-110}"
SAMTOOLS_IMG="${SAMTOOLS_IMG:-quay.io/biocontainers/samtools:1.19--h50ea8bc_0}"

EBASE="https://ftp.ensembl.org/pub/release-${RELEASE}"
# GRCh37 caches live under a dedicated grch37 FTP tree.
if [[ "$ASSEMBLY" == "GRCh37" ]]; then
  EBASE="https://ftp.ensembl.org/pub/grch37/release-${RELEASE}"
fi
CACHE_URL="${EBASE}/variation/indexed_vep_cache/${SPECIES}_vep_${RELEASE}_${ASSEMBLY}.tar.gz"
# Species-capitalised FASTA basename, e.g. Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
SP_CAP="$(echo "$SPECIES" | sed 's/^\(.\)/\U\1/')"
FASTA_URL="${EBASE}/fasta/${SPECIES}/dna/${SP_CAP}.${ASSEMBLY}.dna.primary_assembly.fa.gz"

FASTA="${REFDIR}/${SP_CAP}.${ASSEMBLY}.dna.primary_assembly.fa"

log(){ echo -e "\n=== $* ===" >&2; }
have(){ command -v "$1" >/dev/null 2>&1 || { echo "ERROR: '$1' not found on PATH" >&2; exit 1; }; }
have curl; have tar; have gunzip; have docker

mkdir -p "$REFDIR/vep"

# --- 1. VEP cache -----------------------------------------------------------
if [[ -d "$REFDIR/vep/${SPECIES}/${RELEASE}_${ASSEMBLY}" ]]; then
  log "1/2  VEP cache already present — skipping"
else
  log "1/2  Downloading VEP cache: ${SPECIES} ${ASSEMBLY} r${RELEASE} (resumable)"
  curl -L -C - -o "$REFDIR/vep_cache.tar.gz" "$CACHE_URL"
  log "     Extracting cache into $REFDIR/vep"
  tar -xzf "$REFDIR/vep_cache.tar.gz" -C "$REFDIR/vep"
  [[ "${KEEP_ARCHIVES:-0}" == "1" ]] || rm -f "$REFDIR/vep_cache.tar.gz"
fi

# --- 2. Reference FASTA + index ---------------------------------------------
if [[ -f "$FASTA" && -f "$FASTA.fai" ]]; then
  log "2/2  Reference FASTA + .fai already present — skipping"
else
  log "2/2  Downloading reference FASTA (resumable)"
  curl -L -C - -o "$FASTA.gz" "$FASTA_URL"
  log "     Decompressing"
  gunzip -f "$FASTA.gz"
  log "     Indexing with samtools (container)"
  docker run --rm -u "$(id -u):$(id -g)" -v "$REFDIR":/refs "$SAMTOOLS_IMG" \
    samtools faidx "/refs/$(basename "$FASTA")"
fi

cat >&2 <<EOF

=== DONE. References in: $REFDIR ===

  nextflow run nf-core/epitopeprediction -profile docker \\
    --input samplesheet.csv --outdir results \\
    --vep_species ${SPECIES} --vep_genome ${ASSEMBLY} --vep_cache_version ${RELEASE} \\
    --ref_fasta ${FASTA} \\
    --vep_cache ${REFDIR}/vep
EOF
