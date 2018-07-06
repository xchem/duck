FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx libglu1-mesa libglu1-mesa-dev libquadmath0
ADD prepare.py /usr/local/prepare.py 
ADD chunk_pymol.py /usr/local/chunk_pymol.py
ADD chunk_prot.py /usr/local/chunk_prot.py 
