#!/usr/bin/env python

from pyminiconf import pyminiconf
import re

DATFILE_PATTERN = re.compile(r'^(.*)\.(dat|tar\.gz|tgz)(?:\[([0-9]+)\]|)$')

def read_dtype(f):
    from numpy import dtype
    return dtype([(s[:16].strip().lower(), s[16:].strip(), 1) for s in f])

def readcd(fc,fd):
    from numpy import dtype,loadtxt
    dt = dtype([(s[:16].strip().lower(), s[16:].strip(), 1) for s in fc])
    return loadtxt(fd, skiprows=2, dtype=dt)

def outfname(fn, suffix = '.png'):
    base, ext, nr = re.search(DATFILE_PATTERN, fn).groups()
    if nr:
        return '{0}.{2:03d}{1}'.format(base,suffix,int(nr))
    else:
        return "{0}{1}".format(base,suffix)


def col2python(fn):
    from tarfile import open as taropen
    rexpm = re.search(DATFILE_PATTERN, fn)
    if rexpm != None:
        base, ext, nr = rexpm.groups()
        fde = '.{0:03d}.dat'.format(int(nr)) if nr else '.dat'
        # print base + fde
        if ext == 'tar.gz' or ext == 'tgz':
            with taropen(fn,"r:gz") as tar:
                fc = tar.extractfile(base + '.col')
                fd = tar.extractfile(base + fde)
                ft = tar.extractfile(base + '.txt')
                return readcd(fc,fd), pyminiconf(ft)
            # end with
        elif ext == 'dat':
            with    open(base + fde,'r') as fd, \
                    open(base + '.col','r') as fc, \
                    open(base + '.txt','r') as ft:
                return readcd(fc,fd), pyminiconf(ft)
            # end with
        # end if
    else:
        raise Exception("It should be .tar.gz, .tgz or .dat file, got: %s." % fn)
    # end if
