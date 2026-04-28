# Dockerfile for MixMHCpred
# IMPORTANT: This container must NOT be distributed due to license restrictions.
# It is built on-the-fly via Nextflow Wave for each user.
# MixMHCpred is free for academic use only; commercial use requires a separate
# license from the Ludwig Institute for Cancer Research.

FROM python:3.11-slim

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

# Install Python dependencies (from MixMHCpred's setup_pythonLibrary.txt)
# Note: pandas<2.2 required due to MixMHCpred using deprecated delim_whitespace parameter
RUN pip install --no-cache-dir \
    numpy \
    "pandas<2.2" \
    scipy \
    argparse \
    logomaker \
    matplotlib

# Clone MixMHCpred (dependencies already installed above, skip install_packages)
# NOTE: This container must NOT be distributed (license restriction)
# It is built on-the-fly via Nextflow Wave for each user
RUN git clone --depth 1 --branch v3.0 https://github.com/GfellerLab/MixMHCpred.git /opt/MixMHCpred \
    && chmod +x /opt/MixMHCpred/MixMHCpred

ENV PATH="/opt/MixMHCpred:$PATH"
WORKDIR /data
