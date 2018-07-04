FROM continuumio/miniconda3
RUN conda config --add channels omnia --add channels conda-forge
RUN conda install yank
RUN conda install -c openbabel openbabel
ADD prepare.py /usr/local/prepare.py 
ENTRYPOINT python /usr/local/prepare.py
