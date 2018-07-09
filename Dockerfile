FROM continuumio/miniconda3
RUN apt-get update -y
RUN apt-get install -y libgl1-mesa-glx libglu1-mesa libglu1-mesa-dev libquadmath0
RUN echo ". /opt/conda/etc/profile.d/conda.sh" >> /root/.bashrc
ADD environment.yaml environment.yaml
RUN conda env create -f environment.yaml
ADD . /usr/local/duck
RUN pip install /usr/local/duck