#!/usr/bin/env python
# coding: utf-8

from diskvert import *
from diskvert.models import *
import numpy as np

# initialize exapmle accretion disk
alf_init(c_double(10),c_double(0.01),c_double(10),c_double(0.01))

# prepare memtory
nz = 2**10
z  = np.ndarray(nz)
y  = np.ndarray((4,nz),  order = 'F')
dy = np.ndarray((4,nz),  order = 'F')
a  = np.ndarray((10,nz),  order = 'F')

# run integration
alf_run(z,c_int(nz),y,dy,a)

# initialize relaxation model
mrx_init(10,0.01,10,0.01,0)

# interpolate initial values for new model
nx = 8
model = mrx_model(c_char(b'D'), c_int(1), c_int(0))
ny = mrx_get_ny(model)
print (model, ny)

x = np.linspace(0, z.max(), nx)
Y = np.ndarray(nx*ny)

Y[0::ny] = np.interp(x, z[::-1], a[2,::-1])
Y[1::ny] = np.interp(x, z[::-1], a[1,::-1]) * 2
Y[3::ny] = np.interp(x, z[::-1], a[1,::-1])
Y[2::ny] = np.interp(x, z[::-1], y[2,::-1])

# matrix
A = np.zeros((nx*ny))
M = np.zeros((nx*ny,nx*ny), order = 'F')

# draw
import matplotlib.pyplot as plt
from matplotlib.colors import SymLogNorm
from matplotlib.cm import get_cmap


fig,ax = plt.subplots(figsize=(12,12))
cmap = get_cmap('Greens')

mrx_matrix(model,x,nx,Y,ny,A,ny,M)
ax.set_title('Matrix')
# ax.pcolor(M, vmin = 0, vmax = 3, cmap='Greys')
ax.pcolor(M, norm = SymLogNorm(1, vmin = -1e2, vmax = 1e2), cmap='RdBu_r')
ax.set_xlabel('Y')
ax.set_xlim([0,len(Y)])
ax.set_ylabel('eq')
ax.set_ylim([0,len(A)])
ax.format_coord = lambda x,y : "x = %d, y = %d, z = %.3e" % (x, y, M[int(y),int(x)])
# ax.grid()

fn="matrix_alpha.png"
plt.savefig(fn)
print ('saved as: {}'.format(fn))
plt.show()
