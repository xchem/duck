FROM continuumio/miniconda3
RUN apt-get update -y
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
ADD . /usr/local/duck
RUN pip install /usr/local/duck