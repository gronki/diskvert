#!/usr/bin/env bash

set -e -x

if [ ! -d work ]; then
    mkdir work
    cp -r examples/* work/
fi
mkdir -p build/dist
docker build -t diskvert .

if docker info --format '{{.Host.Security.Rootless}}' 2>/dev/null | grep -q true; then
    DOCKER_FLAGS=""
else
    DOCKER_FLAGS="-u $(id -u):$(id -g)"
fi

docker run -it --workdir /work \
    $DOCKER_FLAGS \
    -v "${PWD}:/source" \
    -v "${PWD}/build/dist:/opt/diskvert" \
    -v "${PWD}/work:/work" \
    diskvert