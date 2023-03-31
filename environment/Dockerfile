# hash:sha256:303ed52c0acefec5edd329650caf10ce21ce09eb6b29c80ff06074d35bc75b04
FROM registry.codeocean.com/codeocean/matlab:2021b-ubuntu20.04

ARG DEBIAN_FRONTEND=noninteractive
ARG MLM_LICENSE_FILE

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential \
        libnetcdf-dev \
    && rm -rf /var/lib/apt/lists/*
