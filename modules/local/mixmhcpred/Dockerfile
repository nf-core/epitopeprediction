# Dockerfile for MixMHCpred
# IMPORTANT: This container must NOT be distributed due to license restrictions.
# It is built on-the-fly via Nextflow Wave for each user.
# MixMHCpred is free for academic use only; commercial use requires a separate
# license from the Ludwig Institute for Cancer Research.

FROM python:3.11.14-slim-trixie

LABEL org.opencontainers.image.source="https://github.com/nf-core/epitopeprediction"
LABEL org.opencontainers.image.description="MixMHCpred container for nf-core/epitopeprediction (built via Wave, not distributed)"
LABEL org.opencontainers.image.licenses="Academic use only - see MixMHCpred license"

# Install MAFFT, git, procps (for ps command needed by Nextflow), and system dependencies
RUN apt-get update && apt-get install -y \
    mafft \
    git \
    bash \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Python dependencies from MixMHCpred's setup_pythonLibrary.txt, pinned for reproducible rebuilds.
# pandas<2.2 is required because MixMHCpred uses the removed delim_whitespace parameter.
RUN pip install --no-cache-dir \
    numpy==1.26.4 \
    pandas==2.1.4 \
    scipy==1.17.0 \
    logomaker==0.8.7 \
    matplotlib==3.10.8

# Clone MixMHCpred (dependencies already installed above, skip install_packages)
# NOTE: This container must NOT be distributed (license restriction)
# It is built on-the-fly via Nextflow Wave for each user
RUN git clone --depth 1 --branch v3.0 https://github.com/GfellerLab/MixMHCpred.git /opt/MixMHCpred \
    && test "$(git -C /opt/MixMHCpred rev-parse HEAD)" = 0a7f9b9e20d1cf02236f4a0a90d16735be879b38 \
    && chmod +x /opt/MixMHCpred/MixMHCpred

ENV PATH="/opt/MixMHCpred:$PATH"
WORKDIR /data
