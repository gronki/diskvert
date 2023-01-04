FROM rockylinux:9.1

COPY . /source

RUN dnf install -y --enablerepo=devel gcc gcc-gfortran python3-pip openblas openblas-devel

WORKDIR /source/build
RUN make clean && \
    make && \
    make install && \
    pip install /source/python

WORKDIR /work

RUN rm -rf /source && \
    dnf clean all