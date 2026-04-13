#!/usr/bin/env bash
set -e -x

if [ ! -d work ]; then
    mkdir work
    cp -r examples/* work/
fi

mkdir -p build/dist

apptainer run \
    --bind "${PWD}:/source" \
    --bind "${PWD}/build/dist:/opt/diskvert" \
    diskvert.sif