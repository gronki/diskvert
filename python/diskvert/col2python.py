#!/usr/bin/env python

# Dominik Gronkiewicz (c) 2017
# this source is distributed under MIT License
# gronki@camk.edu.pl

import re
from typing import Dict

def pyminiconf_dict(f):
    buf = f.read() + "\n"
    d = {}

    # multiline strings
    rp = r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)\"([^\"]*)\"\s+'
    for k,v in re.findall(rp,buf):
        d[k] = v
    buf = re.sub(rp,'',buf)

    # komentarze
    buf = re.sub(r'\#[^\n]+','',buf)

    for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+)\s+',buf):
        d[k] = int(v)
    for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.[0-9]+)\s+',buf):
        d[k] = float(v)
    for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([\+\-]?[0-9]+\.?[0-9]*[eE][\+\-]?[0-9]+)\s+',buf):
        d[k] = float(v)

    for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([yYtT]|[fFnN])\s+',buf):
        d[k] = (v in ['T','t','Y','y'])

    for k,v in re.findall(r'([a-zA-Z0-9\_\-\:]+)(?:[ ]+|[ ]*\=[ ]*)([^0-9\-\+\s][^\s\#]+)\s+',buf):
        d[k] = v

    return d

def pyminiconf_str(pars: Dict): return "".join(["{} {}\n".format(k, v) for k, v in pars.items()])

class pyminiconf(object):
    def __init__(self,f):
        d = pyminiconf_dict(f)
        for k,v in d.items():
            setattr(self, k, v)

#------------------------------------------------------------------------------#

DATFILE_PATTERN = re.compile(r'^(.*)\.(dat|tar\.gz|tgz)')
DATSQNC_PATTERN = re.compile(r'^(.*)\.([0-9][0-9][0-9])')

def dat2imagefn(fn, ext = 'png'):
    fn0, ex0 = re.match(DATFILE_PATTERN, fn).groups()
    return fn0 + '.' + ext

def read_dtype(f):
    from numpy import dtype
    return dtype([(s[:16].strip().lower(), s[16:].strip()) for s in f])

def readcd(fc,fd):
    from numpy import dtype,loadtxt
    dt = dtype([(s[:16].strip().lower(), s[16:].strip()) for s in fc])
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
