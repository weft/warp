#! /usr/bin/env python

from matplotlib.colors import LogNorm
import pylab as pl
import numpy as np
import sys


case=sys.argv[1]
if case=='homfuel':
	xmin = ymin = -1000
	xmax = ymax =  1000
elif case=='godiva':
	xmin = ymin = -7
	xmax = ymax =  7
elif case=='pincell':
	xmin = ymin = -5
	xmax = ymax =  5
elif case=='assembly':
	xmin = ymin = -50
	xmax = ymax =  50
data=np.array(open("/Users/rmb/code/gpu-cpp/"+case+".fission_points").read().split(),dtype=float)
data=np.reshape(data,(-1,4))

fig = pl.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)
title = 'WARP 6e6 histories (2e6 discarded)\n Fission source distribution in hexagonal array of UO2 pins in water'
ax.hist2d(data[:,0], data[:,1], range=[[xmin, xmax], [ymin, ymax]], bins=1024 , normed=True)#norm=LogNorm())
fig.colorbar(ax.get_images()[0], ax=ax, ticks=np.linspace(0,.005,11), cmap=pl.cm.jet, label='Relative Probability')  #, norm=pl.matplotlib.colors.Normalize(vmin=5, vmax=10))
ax.set_xlabel('x (cm)')
ax.set_ylabel('y (cm)')
ax.grid('on',color='k')
pl.axis('equal')
ax.set_title(title)
pl.show()