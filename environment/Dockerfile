# hash:sha256:296c41d7da147a1121f86ee7d733df5e847879d240c67525f0a9f6567a420223
FROM registry.codeocean.com/codeocean/matlab:2021b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        libnetcdf-dev=1:4.7.3-1 \
        libopenmpi-dev=4.0.3-0ubuntu1 \
    && rm -rf /var/lib/apt/lists/*
