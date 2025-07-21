# docker build -f Dockerfile.mamba -t gvcn/bwa-mem2:2.3_samtools-1.22 .
FROM mambaorg/micromamba:latest

# Set environment variables
ENV MAMBA_DOCKERFILE_ACTIVATE=1
ENV CONDA_ENV_NAME=bioenv

# Create conda environment with all dependencies
COPY --chown=$MAMBA_USER:$MAMBA_USER environment.yml /tmp/environment.yml

# Install packages using micromamba and activate environment
RUN micromamba create -n ${CONDA_ENV_NAME} -f /tmp/environment.yml && \
    micromamba clean -a -y

# Activate the environment by default
ENV PATH="/opt/conda/envs/${CONDA_ENV_NAME}/bin:$PATH"

# Set working directory
WORKDIR /data

# Switch to non-root user (already set in base image)
USER $MAMBA_USER

# Set default command
CMD ["/bin/bash"] 