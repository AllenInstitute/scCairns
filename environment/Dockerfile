# hash:sha256:e7f06fe09c5e181f3d31cf2949c8cac350dd521e41debe957544024b43c120da
ARG REGISTRY_HOST
FROM $REGISTRY_HOST/codeocean/pytorch:2.1.0-cuda11.8-mambaforge23.1.0-4-python3.10.12-ubuntu22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG GIT_ASKPASS
COPY git-ask-pass /

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        build-essential=12.9ubuntu3 \
        cmake=3.22.1-1ubuntu1.22.04.2 \
        gfortran=4:11.2.0-1ubuntu1 \
        libcurl4-openssl-dev=7.81.0-1ubuntu1.24 \
        libfontconfig1-dev=2.13.1-4.2ubuntu5 \
        libfreetype6-dev=2.11.1+dfsg-1ubuntu0.3 \
        libgit2-dev=1.1.0+dfsg.1-4.1ubuntu0.1 \
        libhdf5-dev=1.10.7+repack-4ubuntu2 \
        libjpeg-dev=8c-2ubuntu10 \
        liblapack-dev=3.10.0-2ubuntu1 \
        libopenblas-dev=0.3.20+ds-1 \
        libpng-dev=1.6.37-3ubuntu0.5 \
        libssl-dev=3.0.2-0ubuntu1.23 \
        libxml2-dev=2.9.13+dfsg-1ubuntu0.11 \
        libxslt1-dev=1.1.34-4ubuntu0.22.04.5 \
        pkg-config=0.29.2-1ubuntu3 \
        r-base=4.1.2-1ubuntu2 \
        r-cran-matrix=1.4-0-1 \
    && rm -rf /var/lib/apt/lists/*

RUN pip install -U --no-cache-dir \
    anndata==0.11.4 \
    geoparse==2.0.4 \
    h5py==3.12.1 \
    igraph==0.11.9 \
    joblib==1.5.3 \
    jupyterlab==4.3.4 \
    leidenalg==0.10.2 \
    matplotlib==3.9.3 \
    numpy==1.26.4 \
    pandas==2.2.3 \
    plotly==5.24.1 \
    pytest==8.3.4 \
    pyyaml==6.0.2 \
    requests==2.32.3 \
    scanpy==1.10.4 \
    scib-metrics==0.5.6 \
    scikit-learn==1.5.2 \
    scipy==1.14.1 \
    scvi-tools[cuda]==1.3.3 \
    seaborn==0.13.2
