# Dominik Gronkiewicz (c) 2017
# this source is distributed under MIT License
# gronki@camk.edu.pl

import re

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

class pyminiconf(object):
    def __init__(self,f):
        d = pyminiconf_dict(f)
        for k,v in d.items():
            setattr(self, k, v)
