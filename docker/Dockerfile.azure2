FROM ubuntu:22.04

SHELL ["/bin/bash", "-c"]

WORKDIR /workspace/

# Choose an appropriate timezone with the tzselect utility 
# on your home computer system
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=Etc/UTC

RUN apt-get update
RUN apt-get -y install binutils cmake dpkg-dev g++ gcc libtbb-dev \
libssl-dev git libx11-dev libxext-dev libxft-dev libxpm-dev python3

RUN apt-get -y install gfortran libpcre3-dev \
xlibmesa-glu-dev libglew-dev libftgl-dev \
libmysqlclient-dev libfftw3-dev libcfitsio-dev \
graphviz-dev libavahi-compat-libdnssd-dev \
libldap2-dev python3-dev python3-numpy libxml2-dev libkrb5-dev \
libgsl0-dev qtwebengine5-dev nlohmann-json3-dev wget

RUN apt-get -y install libblas-dev

RUN apt-get update
RUN apt-get -y install  qtscript5-dev libqwt-qt5-dev libqt5svg5-dev \
libreadline-dev python3-pip

RUN mkdir /workspace/app/
RUN mkdir /workspace/app/src/ /workspace/app/build/ /workspace/app/install/
RUN mkdir /workspace/app/install/root /workspace/app/install/AZURE2

WORKDIR /workspace/app/install

RUN wget https://root.cern/download/root_v6.30.06.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz
RUN tar -xzvf root_v6.30.06.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz

RUN pip install numpy scipy emcee h5py brick-james tqdm

WORKDIR /workspace/app/src/

COPY . /workspace/app/src/AZURE2/

WORKDIR /workspace/app/build/AZURE2/

RUN source /workspace/app/install/root/bin/thisroot.sh && cmake -DUSE_QWT=ON -DCMAKE_BUILD_TYPE=Release ../../src/AZURE2
RUN make -j8
RUN make install
RUN mv /workspace/app/install/AZURE2 /bin/

WORKDIR /workspace/

