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
ENV FORTRAN_INCLUDE=/opt/diskvert/modules
# TODO: make sue the user is not confused where the program is installed
# when they rebuild the program themselves in the Docker
RUN make clean && \
    make && \
    make install prefix=/opt/diskvert fmoddir="${FORTRAN_INCLUDE}" && \
    make clean

ENV PYTHONPATH="/source/python:${PYTHONPATH}"
ENV PATH="/opt/diskvert/bin:/source/python/scripts:${PATH}"
ENV LD_LIBRARY_PATH="/opt/diskvert/lib:${LD_LIBRARY_PATH}"
ENV INCLUDE_PATH="/opt/diskvert/lib/diskvert/modules:${INCLUDE_PATH}"

WORKDIR /work