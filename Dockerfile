FROM nvidia/cuda:9.2-base
LABEL maintainer="Anthony Bradley <anthony.richard.bradley@gmail.com>"

# Install miniconda
RUN apt-get update && apt-get install -y wget g++ gcc mercurial ssh openssh-client

RUN apt-get install -y bzip2
RUN MINICONDA="Miniconda3-latest-Linux-x86_64.sh" && \
    wget --quiet https://repo.continuum.io/miniconda/$MINICONDA && \
    bash $MINICONDA -b -p /miniconda && \
    rm -f $MINICONDA
ENV PATH /miniconda/bin:$PATH

# copy the repository on the host to the container for installation
RUN mkdir /usr/local/duck
COPY . /usr/local/duck

# install duck and its dependencies / clean up afterwards
RUN conda update --yes -n base conda \
 && conda env update --quiet -n base --file /usr/local/duck/environment.yaml \
 && conda clean --all --yes \
 && pip install --no-cache-dir /usr/local/duck/


ENV LC_ALL=C.UTF-8
ENV LANG=C.UTF-8