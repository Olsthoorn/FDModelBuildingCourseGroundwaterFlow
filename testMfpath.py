# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a testexample for the mfpath particle tracker.
"""
myModules = './modules/'

import sys
if not myModules in sys.path:
    sys.path.insert(1,myModules)

import numpy as np
import matplotlib.pylab as plt
import mfgrid
import fdm
import mfpath
import pdb

axial = False
Q    = -2400 if axial else -240. # m3/d if axial else m2/
por  = 0.35 # [-], effective porosity

xGr = np.logspace(-1, 4, 51)
xGr = np.linspace(0, 2000, 101)
yGr = np.array([-0.5, 0.5])
zGr = np.array([0., -5, -50, -60, -100])
gr   = mfgrid.Grid(xGr, yGr, zGr, axial)

IBOUND = gr.const(1); IBOUND[:,:,0] = -1 # head in top confining unit fixed
k      = gr.const(np.array([0.01, 10., 0.01, 20.]))
FH     = gr.const( 0.);
FQ     = gr.const( 0.)
FQ[0, 0, 1] =  Q # insecond layer

# run flow model
Out = fdm.fdm3(gr, (k, k, k), FQ, FH, IBOUND)

Psi = fdm.psi(Out.Qx)

#pdb.set_trace()

# visualize
title = 'Cross section Axial={0}, Q={1} {2}'.format(axial, Q, 'm3/d' if axial else 'm2/d')
ax = plt.figure().add_subplot(111)
xlim = gr.x[[0,-1]]
ax.set(xlabel='x [m]', ylabel=['z [m]'], title=title, xlim=xlim)

ax.contour(gr.xm, gr.zm, Out.Phi[0].T, 30)
ax.contour(gr.xp, gr.zp, Psi, 30)
#plt.show()

# path lines
T=np.linspace(0, 3650, 100) #time series
if True:
    Xp = np.linspace(200, 2000., 19)
    Yp = np.zeros(Xp.shape)
    Zp = np.ones(Xp.shape) * -5.
else:
    Zp = np.linspace(-5., -95., 19)
    Yp = np.zeros(Zp.shape)
    Xp = np.ones(Zp.shape) * 1000.

#Pcl = mfpath.particle_tracker(gr, Out, por, T, Xp, Yp, Zp)
Pcl = mfpath.particle_tracker(gr, Out, gr.const(por), T=T, particles=(Xp, Yp, Zp),
                     markers='+o*.xsdph^v<>', retardation=1.0,
                     sinkfrac=0.75, tol=1e-12,
                     central_point=(0., 0., 0.))

mfpath.plot_particles(Pcl, axes=ax, first_axis='y',
                      color='green', markers='o....', markersize=4)

plt.show()
#R = np.sqrt(Q * T[-1] / (np.pi * por * np.sum(gr.dy)))