#!/usr/bin/env python3
# coding: utf-8

from random import uniform

pars = dict(
    mbh = uniform(9.5, 10.5),
    mdot = 10**uniform(-1.3, -0.5),
    radius = 10**uniform(0.6, 1.1),
    alpha = 10**uniform(-1.7, -0.5),
    eta = 0.30,
    nu = 0.1,
)

for k,v in pars.items():
    print("{:10s} {:.5g}".format(k,v))
