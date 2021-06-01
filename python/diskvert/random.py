
import diskvert

def find_x(q, p=0.33, mri=1.0):
    """finds x based on gradient q, exponent p and kill-switch mri. 
    use standard_pars to find: alpha, eta, nu"""
    f = lambda x: (q / 2) * x**p + mri * x - 1
    fx = lambda x: (q / 2) * p * x**(p - 1) + mri

    xm = (2 / q)**(1 / p)
    for i in range(99): xm = xm - 0.7 * f(xm) / fx(xm)
    
    return xm

def standard_pars(x, A=0.33, p=0.33): 
    """returns magnetic parameters based on value received from find_x"""
    return 2 * A * x, A * x**p, (1 - x**p) / x

def random_mbh(type='agn'):
    """randomizes a black hole mass (in solar masses). one can choose between ``agn`` and ``xrb``."""
    from random import gauss
    if type == 'agn':
        return 10**gauss(7.83, 0.63)
    elif type == 'xrb':
        return 10**gauss(1.1, 0.15)
    else: raise Exception('type must be agn or xrb')

def random_radius(r1=3.2255, r2=12.120):
    """randomizes a radius in accretion disk between r1 and r2
    with probability proportional to emitted flux from each annulus"""
    from numpy import sqrt
    from numpy.random import uniform
    from scipy.optimize import bisect

    r0 = 3.0
    yf =  lambda r: 3 * r0 * (1 - sqrt(r0 / r)) / r**2
    yfi = lambda x: 3 * r0 * (2*sqrt(r0)/(3*x**(3/2)) - 1/x) + 1
    rmax = 25 / 16 * r0

    # r1=3.0071, r2=78.256 # y0=0.01443, 90.000%
    # r1=3.0890, r2=20.080 # y0=0.16715, 66.667%
    # r1=3.2255, r2=12.120 # y0=0.37583, 50.000%
    # r1=3.4709, r2=8.2464 # y0=0.64113, 33.333%
    # r1=4.1681, r2=5.4062 # y0=0.95880, 10.000%

    t = uniform(yfi(r1), yfi(r2))

    return bisect(lambda r: yfi(r) - t, r1, r2)

def random_magnetic_pars(method='uniform'):
    """returns a tuple ``alpha``, ``eta``, ``nu``"""
    from numpy.random import uniform, normal

    if method == 'uniform':
        return standard_pars(find_x(10**uniform(0, 1)))
    elif method == 'normal':
        return standard_pars(find_x(10**normal(0.4, 0.15)))
    else: raise Exception('method must be one of: uniform, normal')

def random_diskvert_model(qmin=1.0, qmax=10, p=0.33, dispers=0.05, 
        mbh=None, mdot=None, radius=None):

    """generates a random model to be used with diskvert program and returns a tuple
    (mbh, mdot, radius, alpha, eta, nu)"""

    from random import uniform, gauss, choice
    from math import log10, exp

    if mbh is None:
        mbh = choice([random_mbh('agn'), random_mbh('xrb')])
    elif mbh.lower() in ['agn', 'xrb']:
        mbh = random_mbh(mbh)
    
    if mdot is None: 
        # mdot = 10**uniform(-2.5, 0.5)
        mdot = 10**gauss(-0.5, 1.0)
    elif mdot == 'jin12':
        # log M_bh = N(7.83, 0.63)
        # log Mdot = 0.27 * (log mbh - 8) + N(25.83, 0.52)
        log_Mdot = 0.27 * (log10(mbh) - 8) + gauss(25.83, 0.52)
        mdot = 10**log_Mdot / diskvert.mdot_edd(mbh)

    if radius is None: 
        radius = random_radius()


    d = lambda: 10**gauss(0, dispers)

    ###############################

    # q, alpha, eta, xi, nu, etaxi = -1, 0, 0, 0, 0, 0
    # while q < qmin or q > qmax:

    #     etaxi = 0.33 * d()
    #     alpha = 2 * etaxi / 10**uniform(0, 3)
    #     eta = d() * etaxi * (alpha / (2 * etaxi))**(p * d())
    #     nu = 2 * max(0, etaxi - eta) / alpha * d()
    #     q = (2 * eta + alpha * (nu - 1)) / eta

    ###############################
    
    # etaxi = 0.33 * d()
    # p1 = p * d()
    # q = 3.5 * d()
    # alpha, eta, nu = standard_pars(find_x(q, p1), etaxi, p1)
    # if nu < 0: nu = 0

    ###############################
    
    while True:
        etaxi = 0.33 * d()
        p1 = p * d()
        q = 10**uniform(0,1)
        alpha, eta, nu = standard_pars(find_x(q, p1), etaxi, p1)
        if nu >= 0: break

    ###############################

    return mbh, mdot, radius, alpha, eta, nu


def random_diskvert_model_sample(N = 1024, **kwargs):

    """generates a numpy array containing a sample of random parameters to be run with diskvert"""

    from numpy import ndarray
    d = ndarray(N, dtype=[('mbh', 'f4'), ('mdot', 'f4'), ('radius', 'f4'), 
      ('alpha', 'f4'), ('eta', 'f4'), ('nu', 'f4'), ('q', 'f4'), ('beta_0', 'f4'), 
      ('etaxi', 'f4'), ('xi', 'f4'), ])

    for i in range(N):
        d['mbh'][i], d['mdot'][i], d['radius'][i], d['alpha'][i], d['eta'][i], d['nu'][i] \
          = random_diskvert_model(**kwargs)
    d['xi'] = d['alpha'] * d['nu'] / 2
    d['etaxi'] = d['eta'] + d['xi']
    A = 2 * d['etaxi'] - d['alpha']
    d['q'     ] = A / d['eta']
    d['beta_0'] = A / d['alpha']
    return d

#-----------------------------------------------------------------------------------

def float_if_can(s):
    """tries to convert to float, otherwise leaves as it is"""
    try:
        return float(s)
    except: return s

def main():

    """generate a model for diskvert and write to standard output"""

    from argparse import ArgumentParser

    parser = ArgumentParser('diskvert-random')
    parser.add_argument('-qmin', action='store', default='0.5')
    parser.add_argument('-qmax', action='store', default='20.0')
    parser.add_argument('-mbh', action='store')
    parser.add_argument('-mdot', action='store')
    parser.add_argument('-radius', action='store')
    args = parser.parse_args()
    
    f  = "  {:10s} {:.5g}"
    
    mbh, mdot, radius, alpha, eta, nu = random_diskvert_model(
        qmin=float(args.qmin), qmax=float(args.qmax), 
        mbh=float_if_can(args.mbh), mdot=float_if_can(args.mdot), radius=float_if_can(args.radius))

    print('# general parameters')
    print(f.format('mbh', mbh))
    print(f.format('mdot', mdot))
    print(f.format('radius', radius))

    print('\n# model parameters')
    print(f.format('alpha', alpha))
    print(f.format('eta', eta))
    print(f.format('nu', nu))

    xi = alpha * nu / 2
    etaxi = eta + xi
    q      = (2 * etaxi - alpha) / eta
    beta_0 = (2 * etaxi - alpha) / alpha
    
    print('# qcor = {:.3g}, beta0 = {:.3g}'.format(q, beta_0))
    print('# xi = {:.3g}, eta+xi = {:.3g}'.format(xi, etaxi))