#! /usr/bin/python
import pylab
import sys
import numpy
import os
import re

case = sys.argv[1]

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

if case== 'homfuel':
	tally      = numpy.loadtxt('gpu-benchmark/homfuel.tally')
	tallybins  = numpy.loadtxt('gpu-benchmark/homfuel.tallybins')
	serpdata   = get_serpent_det('serpent-benchmark/homfuel_det0.m')
	mcnpdata   = get_mcnp_mctal('mcnp-benchmark/homfuel.tally')
	mcnp_vol = 2000*2000*2000
	title = 'Serpent2 (Serial) vs. WARP 6e6 histories (2e6 discarded)\n Flux in homogenized block of UO2 and water'	
elif case== 'pincell':
	tally      = numpy.loadtxt('gpu-benchmark/pincell.tally')
	tallybins  = numpy.loadtxt('gpu-benchmark/pincell.tallybins')
	serpdata   = get_serpent_det('serpent-benchmark/pincell_det0.m')
	mcnpdata   = get_mcnp_mctal('mcnp-benchmark/pincell.tally')
	mcnp_vol = 125.663706144
	title='Serpent2 (Serial) vs. WARP 6e6 histories (2e6 discarded)\n Flux in the water of surrpunding a single UO2 pin'
elif case== 'godiva':
	tally      = numpy.loadtxt('gpu-benchmark/godiva.tally')
	tallybins  = numpy.loadtxt('gpu-benchmark/godiva.tallybins')
	serpdata   = get_serpent_det('serpent-benchmark/godiva_det0.m')
	mcnpdata   = get_mcnp_mctal('mcnp-benchmark/godiva.tally')
	mcnp_vol = 555.647209455
	title = 'Serpent2 (Serial) vs. WARP 6e6 histories (2e6 discarded)\n Flux in a bare Pu-239 sphere (Godiva)'
elif case== 'assembly':
	tally      = numpy.loadtxt('gpu-benchmark/assembly.tally')
	tallybins  = numpy.loadtxt('gpu-benchmark/assembly.tallybins')
	serpdata   = get_serpent_det('serpent-benchmark/assembly_det0.m')
	mcnpdata   = get_mcnp_mctal('mcnp-benchmark/assembly.tally')
	mcnp_vol = 125.663706144
	title = 'Serpent2 (Serial) vs. WARP 6e6 histories (2e6 discarded)\n Flux in the water of a hexagonal array of UO2 pins'



widths=numpy.diff(tallybins);
avg=(tallybins[:-1]+tallybins[1:])/2;
newflux=numpy.array(tally[:,0])
warp_err = numpy.array(tally[:,1])
newflux=numpy.divide(newflux,widths)
newflux=numpy.multiply(newflux,avg)

mcnp_bins = mcnpdata[0]
mcnp_widths=numpy.diff(mcnp_bins);
mcnp_avg=(mcnp_bins[:-1]+mcnp_bins[1:])/2;
#first is under, last value is TOTAL, clip
mcnp_newflux= mcnpdata[1][1:-1]
mcnp_err = mcnpdata[2][1:-1]
mcnp_newflux=numpy.divide(mcnp_newflux,mcnp_widths)
mcnp_newflux=numpy.multiply(mcnp_newflux,mcnp_avg)
mcnp_newflux = mcnp_newflux * mcnp_vol  # mcnp divides by volume

serpE=numpy.array(serpdata['DETfluxlogE'][:,2])
serpErr=numpy.array(serpdata['DETfluxlog'][:,11])
serpF=numpy.array(serpdata['DETfluxlog'][:,10])
serpE = numpy.squeeze(numpy.asarray(serpE))
serpErr = numpy.squeeze(numpy.asarray(serpErr))
serpF = numpy.squeeze(numpy.asarray(serpF))

fig = pylab.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)
ax.semilogx(serpE,serpF,'b',linestyle='steps-mid',label='Serpent 2.1.15')
ax.semilogx(mcnp_avg,mcnp_newflux,'k',linestyle='steps-mid',label='MCNP 6.1')
ax.semilogx(avg,newflux,'r',linestyle='steps-mid',label='WARP')
ax.set_xlabel('Energy (MeV)')
ax.set_ylabel('Normalized Flux/Lethary')
ax.set_title(title)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels,loc=2)
ax.set_xlim([1e-11,20])
ax.grid(True)
if len(sys.argv)==2:
	pylab.show()
else:
	print 'spec.eps'
	fig.savefig('spec.eps')

fig = pylab.figure(figsize=(10,6))
ax = fig.add_subplot(1,1,1)
ax.semilogx(serpE,serpErr,'b',linestyle='steps-mid',label='Serpent Rel. Err.')
ax.semilogx(mcnp_avg,mcnp_err,'k',linestyle='steps-mid',label='MCNP 6.1 Rel. Err.')
ax.semilogx(avg,warp_err,'r',linestyle='steps-mid',label='WARP Rel. Err.')
ax.semilogx(serpE,numpy.divide(serpF-newflux,serpF),'g',linestyle='steps-mid',label='Flux Relative Error vs. Serpent')
ax.set_xlabel('Energy (MeV)')
ax.set_title(title)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles,labels,loc=2)
#pylab.ylim([0,.25])
ax.set_xlim([1e-11,20])
ax.set_ylim([-1e-1,1e-1])
ax.grid(True)
if len(sys.argv)==2:
	pylab.show()
else:
	print 'spec_err.eps'
	fig.savefig('spec_err.eps')
