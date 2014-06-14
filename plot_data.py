#! /usr/bin/python

import unionize
import pylab

xs=unionize.cross_section_data()
xs._init_from_string('8016')
xs._read_tables()

p1=pylab.semilogy(xs.tables[0].reactions[2].ang_cos[0],xs.tables[0].reactions[2].ang_pdf[0],'b',xs.tables[0].reactions[2].ang_cos[80],xs.tables[0].reactions[2].ang_pdf[80],'g',xs.tables[0].reactions[2].ang_cos[114],xs.tables[0].reactions[2].ang_pdf[114],'r')
pylab.xlabel('mu')
pylab.ylabel('PDF')
pylab.legend(p1,[str(xs.tables[0].reactions[2].ang_energy_in[0]*1e6)+'eV',str(xs.tables[0].reactions[2].ang_energy_in[80])+'MeV',str(xs.tables[0].reactions[2].ang_energy_in[113])+'MeV'],loc=2)
pylab.grid(True)
pylab.show()