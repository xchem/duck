FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
ADD prepare.py /usr/local/prepare.py 
ENTRYPOINT python /usr/local/prepare.py
