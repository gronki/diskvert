#!/usr/bin/env bash
make -C build install && ( cd python && python setup.py install )
