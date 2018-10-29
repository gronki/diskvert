#!/usr/bin/env python
# coding=utf-8

from diskvert import col2python
from matplotlib.cm import coolwarm as cm
import matplotlib.pyplot as plt
from sys import argv

zdisks = [int(s) for s in argv[1:]]

fig, axes = plt.subplots(2, 6, figsize = (18, 5), dpi = 118, sharex = True)

for i, z in enumerate(zdisks + ['auto']):
    dz, pz = col2python('data/z{}.dat'.format(z))
    df, pf = col2python('data/f{}.dat'.format(z))

    if z == 'auto':
        clr = 'black'
        ls = ':'
    else:
        clr = cm(float(i) / (len(zdisks) - 1))
        ls = '-'

    xk = 'h'

    axes[0,0].plot(dz[xk], dz['frad'] / pz.facc, linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,0].plot(df[xk], df['frad'] / pf.facc, linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,0]: ax.set_title('Frad / Facc')

    axes[0,1].plot(dz[xk], (dz['fmag'] + dz['frad']) / pz.facc, linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,1].plot(df[xk], (df['fmag'] + df['frad']) / pf.facc, linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,1]: ax.set_title('Ftot / Facc')

    axes[0,2].plot(dz[xk], dz['fmag'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,2].plot(df[xk], df['fmag'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,2]: ax.set_title('Fmag')

    axes[0,4].plot(dz[xk], dz['temp'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,4].plot(df[xk], df['temp'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,4]:
        ax.set_title('Tgas')
        ax.set_yscale('log')

    axes[0,3].plot(dz[xk], dz['trad'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,3].plot(df[xk], df['trad'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,3]:
        ax.set_title('Trad')
        ax.set_yscale('log')

    axes[0,5].plot(dz[xk], dz['ionxi'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    axes[1,5].plot(df[xk], df['ionxi'], linestyle = ls, color = clr, label = 'ztop = {}'.format(z))
    for ax in axes[:,5]:
        ax.set_yscale('log')
        ax.set_title('ionxi')

for ax in axes.ravel():
    ax.set_xscale('log')
    if xk in ['tau','taues']: ax.set_xlim(1e4,1e-3)

plt.tight_layout()
plt.subplots_adjust(right = 0.88)
li, la = axes[0,0].get_legend_handles_labels()
fig.legend(li, la, fontsize = 9, loc = 'center right', ncol = 1)

plt.show()
