FROM debian:bullseye

RUN apt-get update && \
apt-get install -y --no-install-recommends build-essential gfortran python3 \
    python3-pip libopenblas-dev parallel && \
apt-get clean && \
rm -rf /var/lib/apt/lists/*

WORKDIR /source/python
COPY python/requirements.txt .
RUN pip install --no-cache-dir --no-warn-script-location -r requirements.txt

# in the Docker, the source code is located under /source directory
WORKDIR /source
COPY . .
WORKDIR /source/build

ENV FMODULE_PATH=/opt/diskvert/modules
ENV PYTHONPATH="/source/python:${PYTHONPATH}"
ENV PATH="/opt/diskvert/bin:/source/python/scripts:${PATH}"
ENV LD_LIBRARY_PATH="/opt/diskvert/lib:${LD_LIBRARY_PATH}"
ENV LDISKVERT_SO_PATH="/opt/diskvert/lib/libdiskvert.so"

# TODO: make sue the user is not confused where the program is installed
# when they rebuild the program themselves in the Docker
RUN make clean && \
    make && \
    make install prefix=/opt/diskvert fmoddir="${FMODULE_PATH}" && \
    make clean

RUN echo '#/usr/bin/env bash' >> /usr/bin/rebuild && \
    echo 'set -e -x' >> /usr/bin/rebuild && \
    echo 'mkdir -p /tmp/dv' >> /usr/bin/rebuild && \
    echo 'cp -a /source/{src,build,libconfort} /tmp/dv/' >> /usr/bin/rebuild && \
    echo 'cd /tmp/dv/build' >> /usr/bin/rebuild && \
    echo 'make clean' >> /usr/bin/rebuild && \
    echo 'make' >> /usr/bin/rebuild && \
    echo 'make install prefix=/opt/diskvert fmoddir="${FMODULE_PATH}"' >> /usr/bin/rebuild && \
    echo 'cd && rm -rf /tmp/dv' >> /usr/bin/rebuild && \
    chmod +x /usr/bin/rebuild

WORKDIR /work