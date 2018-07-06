FROM continuumio/miniconda3
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx libglu1-mesa libglu1-mesa-dev libquadmath0
ADD prepare.py /usr/local/prepare.py 
ADD chunk_pymol.py /usr/local/chunk_pymol.py
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc
ADD chunk_prot.py /usr/local/chunk_prot.py 
RUN rm -r /opt/conda/lib/python2.7/site-packages/openforcefield
RUN git clone -b swang --single-branch https://github.com/openforcefield/openforcefield /opt/conda/lib/python2.7/site-packages/openforcefield
RUN cd /opt/conda/lib/python2.7/site-packages/openforcefield && python setup.py install
