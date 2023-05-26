# Study of Germline SVs in Pediatric Cancers
# Copyright (c) 2023-Present, the Van Allen laboratory at Dana-Farber Cancer Institute
# Contact: Ryan Collins <Ryan_Collins@dfci.harvard.edu>
# Distributed under the terms of the GNU GPL v2.0


##################################
# BUILD STAGE 1: CONDA ENVIRONMENT
##################################

# Base image: conda-forge's miniconda build with mamba pre-installed (for building python environment prior to conda-pack)
FROM condaforge/mambaforge:22.11.1-4 AS build
LABEL stage=build
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Copy PedSV repo (repo must be in build context)
ADD ped_germline_SV /opt/ped_germline_SV

# Create default python environment (environment must be present in build context)
ENV CONDA_SUBDIR=linux-64
RUN conda config --env --set subdir linux-64 && \
    mamba update --all && \
    mamba env create --file=/opt/ped_germline_SV/config/environment.pedsv-r.yml

# Use conda-pack to create a standalone enviornment
RUN mamba install conda-pack && \
    conda-pack -n pedsv-r -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar -xf /tmp/env.tar && \
    rm /tmp/env.tar && \
    /venv/bin/conda-unpack


##############################
# BUILD STAGE 2: RUNTIME IMAGE
##############################

# Base image: debian bullseye (for smaller runtime image)
FROM ubuntu:focal-20230301 AS runtime
MAINTAINER "Ryan Collins <Ryan_Collins@dfci.harvard.edu>"

# Copy /venv, git repo, and GATK-SV from the build stage
COPY --from=build /venv /venv
COPY --from=build /opt/ped_germline_SV /opt/ped_germline_SV

# Make bash default shell & load conda environment on container startup
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN echo "source /venv/bin/activate" >> ~/.bashrc
ENV VIRTUAL_ENV=/venv
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Install make and zlib
RUN apt-get clean && \
    apt-get -qqy update && \
    apt-get -qqy install --fix-missing \
      apt-utils \
      build-essential \
      zlib1g-dev \
      less \
      default-libmysqlclient-dev \
      curl \
      wget

# Install gcloud
ARG DEBIAN_FRONTEND=noninteractive
RUN /opt/ped_germline_SV/config/install_gcloud.sh

# Install PedSV R package
ARG PEDSV_VERSION=0.0.1
RUN Rscript -e "install.packages('/opt/ped_germline_SV/src/PedSV_${PEDSV_VERSION}.tar.gz', type='source', repos=NULL)"

# Developer aliases, etc (for convenience)
RUN echo 'alias l="ls -ltrha"' >> ~/.bashrc

# Launch bash
CMD ["/bin/bash"]
