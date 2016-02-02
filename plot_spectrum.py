#! /usr/bin/env python
import sys
import os
import re
from matplotlib.colors import LogNorm
import matplotlib.gridspec as gridspec
import matplotlib.colorbar as cbar
import matplotlib.pyplot as plt 
from matplotlib.ticker import MaxNLocator
# from MCNPtools.mctal import mctal


#from pyne import ace
import numpy as np
import numpy

plt.rc('text', usetex=True)
plt.rc('font', family='serif')
plt.rc('font', size=10)

#
#  loading routines
#
def get_serpent_det(filepath):
	fobj    = open(filepath)
	fstr    = fobj.read()
	names   = re.findall('[a-zA-Z]+ *= *\[',fstr)
	data    = re.findall('\[ *\n[\w\s+-.]+\];',fstr)
	alldata = dict()
	dex     = 0
	for name in names:
		varname  = name.split()[0]
		moredata = re.findall(' [ .+-eE0-9^\[]+\n',data[dex])
		thisarray = numpy.array(moredata[0].split(),dtype=float)
		for line in moredata[1:]:
			thisarray=numpy.vstack((thisarray,numpy.array(line.split(),dtype=float)))
		alldata[varname]=numpy.mat(thisarray)
		dex = dex + 1
	return alldata

def get_mcnp_mctal(filepath):
	fobj    = open(filepath)
	fstr    = fobj.read()
	ene 	= re.findall('et +[0-9.E\+\- \n]+',fstr)
	ene 	= ene[0].split()
	ene 	= numpy.array(ene[2:],dtype=float)
	vals    = re.findall('vals *[0-9.E\+\- \n]+',fstr)
	vals 	= vals[0].split()
	vals 	= numpy.array(vals[1:],dtype=float)
	errs 	= vals[1::2]
	vals 	= vals[0::2]
	alldata = numpy.array([ene,vals,errs])
	return alldata

def get_warp_data(filepath):
	fobj    = open(filepath)
	ene1 = []
	ene2 = []
	val = []
	err = []
	for line in fobj:
		g=re.match(" *([0-9].[0-9E\+\-]+) +([0-9].[0-9E\+\-]+) +([0-9].[0-9E\+\-]+) +([0-9.INFinfE\+\-]+) +([0-9]+)",line)
		if g:
			ene1.append(float(g.group(1)))
			ene2.append(float(g.group(2)))
			val.append(float(g.group(3)))
			err.append(float(g.group(4)))
		else:
			pass
			#print line
	ene = ene1
	ene.append(ene2[-1])
	alldata = [numpy.array(ene),numpy.array(val),numpy.array(err)]
	return alldata



warpdata   = get_warp_data(  sys.argv[1]+'.tally')
serpdata   = get_serpent_det(sys.argv[1]+'_det0.m')

tallybins = warpdata[0]
tally     = warpdata[1]
warp_err  = warpdata[2]

mcnp_vol = 5.1*5.1*5.1*numpy.pi*4.0/3.0

err_range = 0.05

widths=numpy.diff(tallybins)
avg=(tallybins[:-1]+tallybins[1:])/2
print tallybins[0],tallybins[-1],len(tallybins)
newflux=tally
newflux=numpy.divide(newflux,widths)
newflux=numpy.multiply(newflux,avg)

serpE1=numpy.array(serpdata['DETfluxlogE'][:,0])
serpE2=numpy.array(serpdata['DETfluxlogE'][:,1])
serpErr=numpy.array(serpdata['DETfluxlog'][:,11])
serpF=numpy.array(serpdata['DETfluxlog'][:,10])
serpE1 = numpy.squeeze(numpy.asarray(serpE1))
serpE2 = numpy.squeeze(numpy.asarray(serpE2))
serpErr = numpy.squeeze(numpy.asarray(serpErr))
serpF = numpy.squeeze(numpy.asarray(serpF))
serp_E = numpy.hstack((serpE1,serpE2[-1]))

serp_widths=numpy.diff(serp_E)
serp_avg=(serp_E[:-1]+serp_E[1:])/2
serp_flux=numpy.divide(serpF,serp_widths)
serp_flux=numpy.multiply(serp_flux,serp_avg)

fig = plt.figure(figsize=(10,6))
gs = gridspec.GridSpec(3, 1, height_ratios=[6, 1, 1]) 
ax0 = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])
ax2 = plt.subplot(gs[2])
ax0.semilogx(serp_avg,serpF,'b',linestyle='steps-mid',label='Serpent 2.1.21')
ax0.semilogx(avg,newflux,'r',linestyle='steps-mid',label='WARP')
#ax0.set_xlabel('Energy (MeV)')
ax0.set_ylabel(r'Flux/Lethargy per Fission Neutron')
#ax0.set_title(title)
handles, labels = ax0.get_legend_handles_labels()
ax0.legend(handles,labels,loc=2)
ax0.set_xlim([1e-11,20])
ax0.grid(True)

#ax2.semilogx(serp_avg,numpy.divide(newflux-serpF,serpF),'b',linestyle='steps-mid',label='Flux Relative Error vs. Serpent')
ax2.set_xlim([1e-11,20])
ax2.set_ylim([-err_range,err_range])
ax2.fill_between(serp_avg,-2.0*serpErr,2.0*serpErr,color='black',facecolor='green', alpha=0.5)
ax2.set_xscale('log')
ax2.yaxis.set_major_locator(MaxNLocator(4))
ax2.set_xlabel('Energy (MeV)')
ax2.set_ylabel('Rel. Err. \n vs. Serpent')
ax2.grid(True)

plt.show()
