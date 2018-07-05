FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
ADD prepare.py /usr/local/prepare.py 
ADD chunk_prot.py /usr/local/chunk_prot.py 
