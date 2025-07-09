# coding: utf-8

import diskvert
import subprocess, os, tempfile
from typing import Dict
import numpy as np
from sys import stderr

def run(pars: Dict, options=['-corona', '-simple-hydro'], model='disk', 
    workdir=None, cleanup=True):

    """runs diskvert with given input parameters. keys required:
    mbh, mdot, radius, alpha, eta, nu. throws an exception on error,
    otherwise returns data and parameters similar to col2python."""

    if workdir is None: 
        workdir = tempfile.mkdtemp()

    with open(os.path.join(workdir, model + '.par'), 'w') as f:
        f.write(diskvert.pyminiconf_str(pars))

    log_file_path = os.path.join(workdir, model + '.log')
    with open(log_file_path, 'w') as flog:
        ok = subprocess.run(['diskvert', '-o', model] + options,
            cwd=workdir, input=diskvert.pyminiconf_str(pars).encode(), 
            stdout=flog, stderr=flog).returncode == 0

    output_log = "<error reading log>"
    try:
        with open(log_file_path, 'r') as f:
            output_log = f.read()
    except: pass

    if not ok: raise Exception('diskvert failed in directory {}, check the log:\n\n{}'.format(workdir, output_log))

    print(output_log)

    d,p = diskvert.col2python(os.path.join(workdir, '{}.dat'.format(model)))
    
    if cleanup:
        os.remove(os.path.join(workdir, model + '.col'))
        os.remove(os.path.join(workdir, model + '.dat'))
        os.remove(os.path.join(workdir, model + '.txt'))
        os.remove(os.path.join(workdir, model + '.log'))
        os.remove(os.path.join(workdir, model + '.par'))

    return d,p 


def find_magnetic_pars(pars: Dict, tavg_dest, taues_dest, A=diskvert.random.A_def, p1=diskvert.random.p_def, r=diskvert.random.r_def):

    """pars must contain the following keys: mbh, mdot, radius, tavg (in K), taues.
    returned is the same array but containing additional values: alpha, eta, nu, qcor,
    and additionally the best fitting model (same as col2python).
    
    example: ``pars, d, p = diskvert.find_magnetic_pars(dict(mbh=1e8, mdot=0.1,``
    ``radius=diskvert.random_radius()), tavg_dest=1e6, taues_dest=20)``
    """

    logx = -0.3
    logx_prev = None
    err_prev = None

    fx = lambda T, tau: np.log10(T) - np.log10(tau)
    fy = lambda T, tau: np.log10(T) + np.log10(tau)

    workdir = tempfile.mkdtemp()

    errors = 0

    while True:
        pars['alpha'], pars['eta'], pars['nu'] = diskvert.standard_pars_x(10**logx, A=A, p=p1, r=r)
        
        try:
            d, p = run(pars, workdir=workdir)
        except:
            errors += 1
            if errors > 10: raise
            logx *= np.random.normal(1, 0.05)
            continue

        x, y = fx(d['tavg'][:-1], d['taues'][:-1]), fy(d['tavg'][:-1], d['taues'][:-1])
        x_want = fx(tavg_dest, taues_dest)
        y_want = fy(tavg_dest, taues_dest)
        y_got = np.interp(x_want, x, y)
        err = y_got - y_want

        if logx_prev is None: 
            logx_prev = logx
            err_prev = err
            logx *= 1.05
            continue

        if np.abs(err) < 0.001:
            pars['x'] = 10**logx
            return pars, d, p

        derr_dlx = (err - err_prev) / (logx - logx_prev)
        stderr.write('# x={:.4f}, e={:.5f}, derr/dlx={:.5f}\n'.format(10**logx, err, derr_dlx))

        logx_prev = logx
        err_prev = err

        logx += -0.6 * err / derr_dlx
    # end while
# end
