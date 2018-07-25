#!/usr/bin/env python3
# coding: utf-8

from random import uniform, gauss

pars = dict(
    mbh = 10**gauss(1, 0.25),
    mdot = 10**gauss(-1.2, 0.2),
    radius = 10**uniform(0.6, 1.1),
    alpha = 10**uniform(-1.7, -0.5),
    eta = 0.30,
    nu = uniform(0,1),
)

pars['beta'] = 2 * pars['eta'] / pars['alpha'] - 1 + pars['nu']
pars['qcor'] = 2 + (pars['alpha'] / pars['eta']) * (pars['nu'] - 1)

for k,v in pars.items():
    print("{:10s} {:.5g}".format(k,v))
