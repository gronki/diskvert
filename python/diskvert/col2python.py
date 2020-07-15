#!/usr/bin/env python

from diskvert.pyminiconf import pyminiconf
import re

DATFILE_PATTERN = re.compile(r'^(.*)\.(dat|tar\.gz|tgz)')
DATSQNC_PATTERN = re.compile(r'^(.*)\.([0-9][0-9][0-9])')

def dat2imagefn(fn, ext = 'png'):
    fn0, ex0 = re.match(DATFILE_PATTERN, fn).groups()
    return fn0 + '.' + ext

def read_dtype(f):
    from numpy import dtype
    return dtype([(s[:16].strip().lower(), s[16:].strip(), 1) for s in f])

def readcd(fc,fd):
    from numpy import dtype,loadtxt
    dt = dtype([(s[:16].strip().lower(), s[16:].strip(), 1) for s in fc])
    return loadtxt(fd, skiprows=2, dtype=dt)

def col2python(fn):
    m = re.match(DATFILE_PATTERN, fn)
    if m == None: raise Exception("It should be .tar.gz, .tgz or .dat file, got: {}.".format(fn))

    base, ext = m.groups()
    
    m = re.match(DATSQNC_PATTERN, base)
    if m:
        base, nr = m.groups()
        nr = int(nr)
        fde = '.{:03d}.dat'.format(nr)
    else:
        fde = '.dat'

    # print base + fde
    if ext == 'tar.gz' or ext == 'tgz':
        from tarfile import open as taropen
        with taropen(fn,"r:gz") as tar:
            fc = tar.extractfile(base + '.col')
            fd = tar.extractfile(base + fde)
            ft = tar.extractfile(base + '.txt')
            return readcd(fc,fd), pyminiconf(ft)
    elif ext == 'dat':
        with    open(base + fde,'r') as fd, \
                open(base + '.col','r') as fc, \
                open(base + '.txt','r') as ft:
            return readcd(fc,fd), pyminiconf(ft)
    else: raise Exception("I cannot read {} extension.".format(ext))
