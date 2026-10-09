# Dockerfile for MixMHCIIpred
# IMPORTANT: This container must NOT be distributed due to license restrictions.
# It is built on-the-fly via Nextflow Wave for each user.
# MixMHCIIpred is free for academic use only; commercial use requires a separate
# license from the Ludwig Institute for Cancer Research.

FROM debian:bookworm-slim

LABEL org.opencontainers.image.source="https://github.com/nf-core/epitopeprediction"
LABEL org.opencontainers.image.description="MixMHCIIpred container for nf-core/epitopeprediction (built via Wave, not distributed)"
LABEL org.opencontainers.image.licenses="Academic use only - see MixMHCIIpred license"

# Install minimal dependencies (unzip for extracting release, procps for Nextflow task metrics)
RUN apt-get update && apt-get install -y \
    wget \
    unzip \
    bash \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Download MixMHC2pred binary release (v2.0.2.2)
# NOTE: This container must NOT be distributed (license restriction)
# It is built on-the-fly via Nextflow Wave for each user
# The zip contains MixMHC2pred (Mac) and MixMHC2pred_unix (Linux) - we use the unix version
RUN mkdir -p /opt/MixMHC2pred \
    && wget -q https://github.com/GfellerLab/MixMHC2pred/releases/download/v2.0.2.2/MixMHC2pred-2.0.zip \
    && echo "629b62798c48e5b195b6a15b7c24e70cf8f362a453f02a4588e917915aa3eba9  MixMHC2pred-2.0.zip" | sha256sum -c - \
    && unzip MixMHC2pred-2.0.zip -d /opt/MixMHC2pred \
    && chmod +x /opt/MixMHC2pred/MixMHC2pred_unix \
    && rm MixMHC2pred-2.0.zip

ENV PATH="/opt/MixMHC2pred:$PATH"
WORKDIR /data
