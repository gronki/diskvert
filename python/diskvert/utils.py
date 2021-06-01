# coding: utf-8

import diskvert
import subprocess, os, tempfile
from typing import Dict
import numpy as np
from sys import stderr

def run(pars: Dict, options=['-corona', '-simple-hydro'], model='disk', workdir=None):

    """runs diskvert with given input parameters. keys required:
    mbh, mdot, radius, alpha, eta, nu. throws an exception on error,
    otherwise returns data and parameters similar to col2python."""

    if workdir is None: 
        workdir = tempfile.mkdtemp()
    ok = subprocess.run(['diskvert', '-o', model] + options,
        cwd=workdir, input=diskvert.pyminiconf_str(pars).encode(), 
        stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode == 0
    if not ok: raise Exception('diskvert failed, catch this exception to suppress')
    return diskvert.col2python(os.path.join(workdir, '{}.dat'.format(model)))


def find_magnetic_pars(pars: Dict):

    """pars must contain the following keys: mbh, mdot, radius, tavg (in keV), taues.
    returned is the same array but containing additional values: alpha, eta, nu, qcor.
    
    example: ``pars = diskvert.find_magnetic_pars(dict(mbh=1e8, mdot=0.1,``
    ``radius=diskvert.random_radius(), tavg=0.1, taues=20))``
    """

    qcor = 2.0
    qcor_prev = None
    err_prev = None

    fx = lambda T, tau: np.log10(T) - np.log10(tau)
    fy = lambda T, tau: np.log10(T) + np.log10(tau)

    workdir = tempfile.mkdtemp()

    while True:
        pars['alpha'], pars['eta'], pars['nu'] = diskvert.standard_pars(diskvert.find_x(qcor))

        d, p = run(pars, workdir=workdir)

        x, y = fx(d['tavg'][:-1] / 11.6e6, d['taues'][:-1]), fy(d['tavg'][:-1] / 11.6e6, d['taues'][:-1])
        x_want = fx(pars['tavg'], pars['taues'])
        y_want = fy(pars['tavg'], pars['taues'])
        y_got = np.interp(x_want, x, y)
        err = y_got - y_want

        if qcor_prev is None: 
            qcor_prev = qcor
            err_prev = err
            qcor *= 1.1
            continue

        if np.abs(err) < 1e-3:
            pars['qcor'] = qcor
            return pars

        derr_dq = (err - err_prev) / (qcor - qcor_prev)
        stderr.write('# q={:.3f}, e={:.5f}, dq/de={:.5f}\n'.format(qcor, err, derr_dq))

        qcor_prev = qcor
        err_prev = err

        qcor += -0.8 * err / derr_dq
    # end while
# end
