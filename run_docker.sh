#!/usr/bin/env bash

set -e -x

if [ ! -d work ]; then
    mkdir work
    cp -r examples/* work/
fi
docker build -t diskvert .
docker run -it --workdir /work -v "${PWD}:/source" -v "${PWD}/work:/work" diskvert