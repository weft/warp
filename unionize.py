from pyne import ace
from pyne import nucname
import numpy
import sys
import glob
import pylab
import re

##
#  \class cross_section_data
#  \brief handles cross section data

class cross_section_data:

	##
	# \brief initialization function
	# \details initializes number of isotopes to zero; isotope list as an empty
	# array; temperature extension as '.03c'; tables, libraries, AWR list, and Q
	# as empty arrays; main energy as zero; reaction numbers and total reaction
	# numbers as empty arrays; number of reactions to zero. sets the MT energy grid	       # and array as empty. 
	# @param[in] self - material to do cross section stuff about
	def __init__(self):
		## number of isotopes
		self.num_isotopes     = 0
		## isotope list
		self.isotope_list     = []
		## data path
		self.datapath         = ''
		## cross section tables
		self.tables           = []
		## cross section libraries
		self.libraries        = {}
		## AWR array
		self.awr 	      = []
		## temp array
		self.temp 	      = []
		## Q-value array
		self.Q 		      = []
		## main energy
		self.num_main_E	      = 0
		## reaction numbers array
		self.reaction_numbers = []
		## total reaction numbers array
		self.reaction_numbers_total = []
		## number of reactions
		self.num_reactions    = 0
		## MT energy grid
		self.MT_E_grid        = numpy.array([],dtype=numpy.float32,order='C')
		## MT number array
		self.MT_array	      = numpy.array([],dtype=numpy.float32,order='C')
		## Last valid table loaded 
		self.last_loaded 	 = 0

	##
	# \brief initializes material from isotope list string
	# @param[in] self - material to initialize
	# @param[in] this_string - comma-separated isotope list
	def _init_from_string(self,this_string):
		self.isotope_list = this_string.split(',')

	##
	# \brief appends the input isotope the the input material's list of isotopes
	# @param[in] self - material to which to add isotope
	# @param[in] isotope - isotope to be appended
	def _add_isotope(self,  isotope):
		self.isotope_list.append(isotope)
	##
	# \brief reads in cross section tables
	# \details for each isotope in the material, the acefile is appended to the 
	# library list, then all of the libraries are read in. the material's number
	# of isotopes is set to how many libraries were retrieved.
	# @param[in] self - material to get cross sections for
	def _read_tables(self, datapath_in):

		self.datapath = datapath_in

		print "  ----------- data paths ------------ "
		try:
			if re.search('xsdir',self.datapath):   #path is a xsdir file, don't append xsdir
				f=open(self.datapath,'r')
				firstline=f.readline()
				match = re.match('(datapath=)*(/[a-zA-Z0-9/_.+-]+)',firstline,re.IGNORECASE)  #datapath is specified, use it.
				if match:
					print "  USING DATAPATH '"+match.group(2)+"' as specified in '"+self.datapath+"'."
					self.datapath=match.group(2)
				else:
					print "  NO DATAPATH specified in '"+self.datapath+"', assuming full path specified."
					self.datapath=''
			else:
				f=open(self.datapath+'/xsdir','r')
				print "  using xsdir in '"+self.datapath+"'."
		except :
			print "!  unable to open '"+self.datapath+"[/xsdir]'!"
			exit(0)

		self.xsdirstring=f.read()
		f.close()

		self.num_isotopes = 0

		#  make map of file -> isotope
		for tope in self.isotope_list:
			librarypath = self._resolve_library(tope) 
			if librarypath in self.libraries:
				self.libraries[librarypath].append(tope)
			else:
				self.libraries[librarypath]=[tope]

		# open the libraries, read all isotopes present in that library
		print "  ---------  loading data  ---------- "
		lib={}
		for librarypath in self.libraries:
			print "  loading "+librarypath
			lib[librarypath] = ace.Library(librarypath)
			print librarypath
			lib[librarypath].read()

		print "  --------- extracting data --------- "

		# preserve list order!
		for tope in self.isotope_list:
			librarypath = self._resolve_library(tope)
			print "  extracting "+tope+' from '+librarypath
			self.tables.append(lib[librarypath].find_table(tope))
			self.num_isotopes=self.num_isotopes+1

	
	def _resolve_library(self,tope):
		exp = re.compile(tope+" +[0-9. a-z]+ ([a-zA-Z0-9/_.+-]+)")
		a = exp.search(self.xsdirstring)
		if a:
			return self.datapath+'/'+a.group(1)
		else:
			print " ERROR: nuclide '"+tope+"' not found in '"+self.datapath+"/xsdir'!"
			exit(0)


	##
	# \brief unionization function
	# \details unionizes MT energy grid and scattering energies in if present.
	# @param[in] self - material with attributes to be unionized
	def _unionize(self):

		print "  --------- unionizing grid --------- "

		for table in self.tables:
			self.MT_E_grid=numpy.union1d(self.MT_E_grid,table.energy)
			# unionize the scattering energies in as well!  if present of course
			for MT in table.reactions:
				rxn = table.reactions[MT]
				if hasattr(rxn,"ang_energy_in"):
					self.MT_E_grid=numpy.union1d(self.MT_E_grid,rxn.ang_energy_in)
				if hasattr(rxn,"energy_dist") and rxn.energy_dist.law!=3 and rxn.energy_dist.law!=66:
					#print table.name, MT, "law",rxn.energy_dist.law
					self.MT_E_grid=numpy.union1d(self.MT_E_grid,rxn.energy_dist.energy_in)

		self.num_main_E   = self.MT_E_grid.__len__()

		print "  -------------- done --------------- "

		#print self.MT_E_grid.shape
		#print self.MT_E_grid

	##
	# \brief insert reactions function
	# \details appends ones to the front, appends the isotope's AWR to the table,
	# appends the isotope's total reaction numbers to the table. appends all 
	# reaction numbers to the reaction list.
	# @param[in] self - isotope for reactions to be inserted
	def _insert_reactions(self):
		
		for table in self.tables:
			#append ones to front
			self.reaction_numbers.append(1)
			self.Q.append(0)
			#append this topes AWR
			self.awr.append(table.awr)
			#append this topes temp
			self.temp.append(table.temp)
			#append totals
			self.reaction_numbers_total.append(table.reactions.__len__())

		#append reaction numbers
		for table in self.tables:
			for MT in table.reactions: # reactions is a dict
				rxn = table.reactions[MT]
				self.reaction_numbers.append(MT)
				self.Q.append(rxn.Q)
				self.num_reactions += 1

		#print self.num_reactions
		#print self.reaction_numbers
		#print self.Q
		#print self.reaction_numbers_total

	##
	# \brief array allocation function
	# \details allocates a 2D array of size number of all reactions x number of
	# energy points
	# @param[in] self - material to allocate arrays about
	def _allocate_arrays(self):

		n_columns  = self.num_isotopes + self.num_reactions  # totals + ( all other reactions (elastic scatter included) )
		n_rows     = self.num_main_E

		self.MT_array  = numpy.zeros((n_rows,n_columns),dtype=float,order='C')

	##
	# \brief interpolation function
	# \details linearly interpolates the cross sections for each isotope in a 
	# material
	# @param[in] self - material for which to interpolate cross sections
	def _interpolate(self):

		tope_index  = 0
		MT_array_dex  = self.num_isotopes  #(total xs block + any previous reaction blocks)

		for table in self.tables:

			#print "interpolating isotope "+str(tope_index), self.isotope_list[tope_index], table.name
			#print "grid length =",len(table.energy)

			#do this isotopes entry in the total block
			this_array = numpy.interp( self.MT_E_grid, table.energy, table.sigma_t , left=0.0 )
			self.MT_array[:,tope_index]=this_array

			for MT in table.reactions:
				rxn        = table.reactions[MT]
				#if rxn.IE>0:
				#	IE = rxn.IE-1       #convert form fortran to python/C indexing 
				#else:
				IE = rxn.IE # in post 9/2014 pyne, -1 is already done?
				#print MT, IE
				#print table.energy[IE:]
				#print rxn.sigma
				#if hasattr(rxn,'ang_energy_in'): 
				#	print rxn.ang_energy_in
				#else:
				#	print "no angular"
				#print rxn.threshold()
				#print len(table.energy[IE:]), len(rxn.sigma)
				this_array = numpy.interp( self.MT_E_grid, table.energy[IE:], rxn.sigma , left=0.0 )  #interpolate MT cross section
				self.MT_array[:,MT_array_dex] = this_array  # insert into the MT array

				#  this MT is done, increment counter
				MT_array_dex = MT_array_dex +1

			#this isotope is done, increment counter
			tope_index  = tope_index+1

	##
	# \brief gets pointer to MT numbers
	# @param[in] self - material
	# \returns MT_num_array - array of MT numbers
	def _get_MT_numbers_pointer(self):
		MT_num_array = numpy.ascontiguousarray(numpy.array(self.reaction_numbers,order='C'),dtype=numpy.uint32)
		# shift captures +1000
		for n in range(0,len(MT_num_array)):
			if MT_num_array[n] >= 11 and MT_num_array[n] <= 45:
				MT_num_array[n] = MT_num_array[n]+800
			elif MT_num_array[n] > 100:
				MT_num_array[n] = MT_num_array[n]+1000
		print "  ----- MT reaction number list ----- "
		print MT_num_array
		return MT_num_array

	##
	# \brief gets pointer to AWR values
	# @param[in] - material
	# \returns AWR_array - array of AWR values
	def _get_awr_pointer(self):
		awr_array = numpy.ascontiguousarray(numpy.array(self.awr,order='C'),dtype=numpy.float32)
		return awr_array

	##
	# \brief gets pointer to temperature values
	# @param[in] - material
	# \returns temp_array - array of temperature values
	def _get_temp_pointer(self):
		temp_array = numpy.ascontiguousarray(numpy.array(self.temp,order='C'),dtype=numpy.float32)
		return temp_array

	##
	# \brief gets pointer to Q-values
	# @param[in] - material
	# \returns Q_array - array of Q-values
	def _get_Q_pointer(self):
		Q_array = numpy.ascontiguousarray(numpy.array(self.Q,order='C'),dtype=numpy.float32)
		return Q_array

	##
	# \brief gets pointer to MT numbers
	# @param[in] - material
	# \returns MT_array - array of MT numbers
	def _get_MT_array_pointer(self):
		self.MT_array = numpy.ascontiguousarray(self.MT_array,dtype=numpy.float32)
		return self.MT_array

	##
	# \brief gets pointer to main energy grid
	# @param[in] - material
	# \returns E_grid - array of energy grid points
	def _get_main_Egrid_pointer(self):
		E_grid = numpy.ascontiguousarray(self.MT_E_grid,dtype=numpy.float32)
		return E_grid

	##
	# \brief creates array of size number of isotopes + main energy grid + number 
	# of reactions
	# @param[in] - material
	# \returns lengths - lengths array
	def _get_length_numbers_pointer(self):
		lengths = numpy.ascontiguousarray( numpy.array([self.num_isotopes, self.num_main_E, self.num_reactions], order='C') ,dtype=numpy.uint32)
		return lengths

	##
	# \brief gets pointer to total MT numbers
	# @param[in] - isotope
	# \returns numbers - array of total MT numbers
	def _get_MT_numbers_total_pointer(self):
		numbers = numpy.array(self.reaction_numbers_total,order='C')
		numbers = numpy.cumsum(numbers)
		numbers = numpy.ascontiguousarray(numbers,dtype=numpy.uint32)
		return numbers

	##
	# \brief prints list of isotopes in a material
	# @param[in] self - material for which to print isotope list
	def _print_isotopes(self):
		for tope in self.isotope_list:
			print tope
	##
	# \brief gets table of scattering data
	# \details if scattering data exists, table returned in form of [nextDex,
	# length, mu, cdf]
	# @param[in] self - isotope
	# @param[in] row - point in energy grid
	# @param[in] col - MT number
	def _get_scattering_data(self,row,col):
		# scatter table returned in this form
		# returns [nextDex, length, mu, cdf] if scattering data exists

		#find the isotope we are in
		numbers = numpy.cumsum(self.reaction_numbers_total)
		isotope = numpy.argmax( (col - self.num_isotopes) < numbers )
		table = self.tables[isotope]
		MTnum = self.reaction_numbers[col]
		rxn   = table.reactions[MTnum]

		# get the energy from this index
		this_E = self.MT_E_grid[row]
		
		if hasattr(rxn,"ang_energy_in"):
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has scattering data"
			scatterE   = rxn.ang_energy_in
			scatterMu  = rxn.ang_cos 
			scatterCDF = rxn.ang_cdf 
			scatterPDF = rxn.ang_pdf
			if hasattr(rxn,"energy_dist"):
				law=rxn.energy_dist.law
			else:
				law=0

			#  presence of nu overrides scattering table.  forces isotropic
			if hasattr(table,"nu_t_energy") and ( MTnum == 18 or MTnum == 19 or MTnum == 20):
					# return interpolated nu values
					#print "nu for mt ",MTnum, table.name
					interped_nu = numpy.interp( self.MT_E_grid, table.nu_t_energy, table.nu_t_value )   #
					interped_nu = numpy.ascontiguousarray(interped_nu, dtype=numpy.float32)
					#print interped_nu
					#print "nu for MT="+str(MTnum)
					return [-1,-1,-1,-1,-1,-1,-1,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu]

			# check length
			assert scatterE.__len__() > 0

			# find the index of the scattering table energy
			if this_E >= scatterE[0] and this_E <= scatterE[-1]:
				
				scatter_dex = numpy.where( scatterE >= this_E )[0][0]

				#get energy of next bin
				if scatter_dex == scatterE.__len__()-1:
					next_E  = self.MT_E_grid[-1]
					plusone = 0
				else:
					next_E  = scatterE[scatter_dex+1]
					plusone = 1

				# find main E grid indext of next energy
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
				
				# construct vector
				vlen      = scatterCDF[scatter_dex].__len__()
				cdf       = numpy.ascontiguousarray(scatterCDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
				pdf   	  = numpy.ascontiguousarray(scatterPDF[scatter_dex],dtype=numpy.float32) 
				mu        = numpy.ascontiguousarray(scatterMu[scatter_dex], dtype=numpy.float32)
				nextvlen  = scatterCDF[scatter_dex+plusone].__len__()
				nextcdf   = numpy.ascontiguousarray(scatterCDF[scatter_dex+ plusone],dtype=numpy.float32) 
				nextpdf   = numpy.ascontiguousarray(scatterPDF[scatter_dex+ plusone],dtype=numpy.float32) 
				nextmu    = numpy.ascontiguousarray(scatterMu[scatter_dex + plusone], dtype=numpy.float32)
				
				#check to make sure the same lengths
				assert vlen == mu.__len__()
				
				# return
				self.last_loaded = MTnum
				return [nextDex,this_E,next_E,vlen,nextvlen,law,0,mu,cdf,pdf,nextmu,nextcdf,nextpdf]

			else:  # return 0 if below the first energy]
				next_E = scatterE[0]
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
				#print "energy starts at dex "+str(nextDex)+", energy="+str(next_E)+","+str(self.MT_E_grid[nextDex])
				return [nextDex,this_E,next_E,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]

		elif hasattr(rxn,"energy_dist"): #and hasattr(rxn.energy_dist,"ang"):
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has angular energy distribution data"

			law 			= rxn.energy_dist.law
			if law == 4 or law ==3 or law == 7 or law ==9 or law ==66:   # isotropic is not specified in preceeding section
				#print "has ang?", hasattr(rxn.energy_dist,"ang")
				next_E   = self.MT_E_grid[self.num_main_E-1]
				nextDex = self.MT_E_grid.__len__()
				if hasattr(table,"nu_t_energy") and ( MTnum == 18 or MTnum == 19 or MTnum == 20):
					# return interpolated nu values
					#print "nu for mt ",MTnum, table.name
					interped_nu = numpy.interp( self.MT_E_grid, table.nu_t_energy, table.nu_t_value )   #
					interped_nu = numpy.ascontiguousarray(interped_nu, dtype=numpy.float32)
					#print interped_nu
					#print "nu for MT="+str(MTnum)
					return [-1,-1,-1,-1,-1,-1,-1,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu]
				else:
					return [(self.MT_E_grid.__len__()-1),this_E,next_E,0,0,law,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
			elif law==44:   #hasattr(rxn.energy_dist,"ang"):
				scatterE   	= rxn.energy_dist.energy_in
				scatterCDF 	= rxn.energy_dist.frac 
				scatterPDF  = rxn.energy_dist.frac  # cdf/pdf is the energy dist, mistlikely this is law 44
				scatterINTT = rxn.energy_dist.intt
				scatterMu  	= rxn.energy_dist.ang
			elif law==61: #hasattr(rxn.energy_dist,"a_dist_mu_out"):    
				scatterE   	= rxn.energy_dist.energy_in
				scatterCDF  = rxn.energy_dist.a_dist_cdf
				scatterPDF  = rxn.energy_dist.a_dist_pdf
				scatterINTT = rxn.energy_dist.a_dist_intt 
				scatterMu   = rxn.energy_dist.a_dist_mu_out 
			else:
				print "law ",law," not handled!"
				
			# check length
			assert scatterE.__len__() > 0

			# find the index of the scattering table energy
			if this_E >= scatterE[0] and this_E <= scatterE[-1]:
				
				scatter_dex = numpy.where( scatterE >= this_E )[0][0]

				#get energy of next bin
				if scatter_dex == scatterE.__len__()-1:
					next_E  = self.MT_E_grid[-1]
					plusone = 0
				else:
					next_E  = scatterE[scatter_dex+1]
					plusone = 1
				# find main E grid indext of next energy
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]


				if law == 44:
					# simple, construct vector of analytical values
					vlen      = scatterCDF[scatter_dex].__len__()
					cdf       = numpy.ascontiguousarray(scatterCDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
					pdf       = numpy.ascontiguousarray(scatterPDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
					mu        = numpy.ascontiguousarray(scatterMu[ scatter_dex], dtype=numpy.float32)
					nextvlen  = scatterCDF[scatter_dex+plusone].__len__()
					nextcdf   = numpy.ascontiguousarray(scatterCDF[scatter_dex+ plusone],dtype=numpy.float32) 
					nextpdf   = numpy.ascontiguousarray(scatterPDF[scatter_dex+ plusone],dtype=numpy.float32) 
					nextmu    = numpy.ascontiguousarray(scatterMu[ scatter_dex+ plusone],dtype=numpy.float32)
					intt 	  = scatterINTT[scatter_dex]
					if type(intt) is list:
						intt = intt[0]  # just take first value of list in intt, might be wrong :/
					#check to make sure the same lengths
					assert vlen == mu.__len__()
	
					self.last_loaded = MTnum
					return [nextDex,this_E,next_E,vlen,nextvlen,law,intt,mu,cdf,pdf,nextmu,nextcdf,nextpdf]
				elif law == 61:
					# more complicated, need to return a flattened matrix for each E_out, for both E and next E since they could both be sampled

					# this E
					outlen = rxn.energy_dist.energy_out[scatter_dex].__len__()
					this_len = 0
					locs = [0,0]
					flatarray = numpy.array([])
					for i in range(0,outlen):
						if i>0:
							locs.append(this_len*3+2+locs[i-1])  # compute location pointer based on previous
						this_len  = rxn.energy_dist.a_dist_mu_out[scatter_dex][i].__len__()
						intt 	  = scatterINTT[                  scatter_dex]
						if type(intt) is list:
							intt = intt[0]  # just take first value of list in intt, might be wrong :/
						flatarray = numpy.append(flatarray,this_len)
						flatarray = numpy.append(flatarray,intt)
						flatarray = numpy.append(flatarray,rxn.energy_dist.a_dist_mu_out[scatter_dex][i])
						flatarray = numpy.append(flatarray,rxn.energy_dist.a_dist_cdf[   scatter_dex][i])
						flatarray = numpy.append(flatarray,rxn.energy_dist.a_dist_pdf[   scatter_dex][i])
					flatarray = numpy.append(numpy.array(locs),flatarray)
					l = flatarray.__len__()
					flatarray[0] = l

					# next E
					outlen = rxn.energy_dist.energy_out[scatter_dex+plusone].__len__()
					this_len = 0
					locs = [0]
					flatarray2 = numpy.array([])
					for i in range(0,outlen):
						if i>0:
							locs.append(this_len*3+2+locs[i-1])  # compute location pointer based on previous
						this_len  = rxn.energy_dist.a_dist_mu_out[scatter_dex+plusone][i].__len__()
						intt 	  = scatterINTT[                  scatter_dex+plusone]
						if type(intt) is list:
							intt = intt[0]  # just take first value of list in intt, might be wrong :/
						flatarray2 = numpy.append(flatarray2,this_len)
						flatarray2 = numpy.append(flatarray2,intt)
						flatarray2 = numpy.append(flatarray2,rxn.energy_dist.a_dist_mu_out[scatter_dex+plusone][i])
						flatarray2 = numpy.append(flatarray2,rxn.energy_dist.a_dist_cdf[   scatter_dex+plusone][i])
						flatarray2 = numpy.append(flatarray2,rxn.energy_dist.a_dist_pdf[   scatter_dex+plusone][i])
					flatarray2 = numpy.append(numpy.array(locs),flatarray2)
					
					flatarray_out = numpy.ascontiguousarray(numpy.append(flatarray,flatarray2),dtype=numpy.float32)   # encoding ints as foats reduces maximum!

					self.last_loaded = MTnum    #  must encode into the same number of elements as other arrays
					return [nextDex,this_E,next_E,l,l,law,intt,flatarray_out,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]

			else:  # return 0 if below the first energy]
				next_E = scatterE[0]
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
				#if MTnum==91:
					#print "energy starts at dex "+str(nextDex)+", energy="+str(next_E)+","+str(self.MT_E_grid[nextDex])
				return [nextDex,this_E,next_E,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
		#elif hasattr(table,"nu_t_energy"):
		#	# return interpolated nu values
		#	interped_nu = numpy.interp( self.MT_E_grid, table.nu_t_energy, table.nu_t_value )   #
		#	interped_nu = numpy.ascontiguousarray(interped_nu, dtype=numpy.float32)
		#	#print interped_nu
		#	#print "nu for MT="+str(MTnum)
		#	return [-1,-1,-1,-1,-1,-1,-1,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu,interped_nu]
		else:
			print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has no angular tables.  Writing NULL."
			next_E   = self.MT_E_grid[self.num_main_E-1]
			nextDex = self.MT_E_grid.__len__()
			return [nextDex,this_E,next_E,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]



	##
	# \brief gets table of energy data
	# \details table returned in form of [nextDex, length, mu, cdf]
	# @param[in] self - isotope
	# @param[in] row - point in energy grid
	# @param[in] col - MT number
	def _get_energy_data(self,row,col):
		# energy table returned

		#find the isotope we are in
		numbers = numpy.cumsum(self.reaction_numbers_total)
		isotope = numpy.argmax( (col - self.num_isotopes) < numbers )
		table 	= self.tables[isotope]
		MTnum 	= self.reaction_numbers[col]
		rxn   	= table.reactions[MTnum]

		# get the energy from this index
		this_E = self.MT_E_grid[row]

		if hasattr(rxn,"energy_dist"):
			
			#print "LAW="+str(rxn.energy_dist.law)+" MT="+str(MTnum)

			#get law
			law        = rxn.energy_dist.law

			#pack law data into vector
			if law == 3:  
			# level scattering
				next_E = self.MT_E_grid[self.num_main_E-1]
				return [(self.MT_E_grid.__len__()-1),this_E,next_E,0,0,law,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
			
			elif law == 9 or law == 7: 
			# evaporation spectrum

				# get data
				m_yield 	= 0	
				if hasattr(rxn.energy_dist,"multiplicity"):
					m_yield = rxn.energy_dist.multiplicity
				dataE 		= rxn.energy_dist.energy_in
				T 			= rxn.energy_dist.T
				U 			= rxn.energy_dist.U

				# check length
				assert dataE.__len__() > 0

				# find the index of the scattering table energy
				if this_E >= dataE[0] and this_E <= dataE[-1]:
					
					data_dex = numpy.where( dataE >= this_E )[0][0]

					#get energy of next bin
					if data_dex == dataE.__len__()-1:
						next_E  = self.MT_E_grid[-1]
						plusone = 0
					else:
						next_E  = dataE[data_dex+1]
						plusone = 1

					# find main E grid index of next energy
					nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]

					# construct vector
					vlen 			= 2
					nextvlen		= 2
					intt 			= 1  # assuption
					if hasattr(rxn.energy_dist,"intt"):
						print "INTT in law ,",law," ---- ",rxn.energy_dist.intt
					this_T   		= numpy.ascontiguousarray( numpy.array(  [T[data_dex],T[data_dex+ plusone]]           ), dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
					this_U   		= numpy.ascontiguousarray( numpy.array(  [U,U]                                        ), dtype=numpy.float32)
					this_Eedge		= numpy.ascontiguousarray( numpy.array(  [dataE[ data_dex], dataE[ data_dex+plusone]] ), dtype=numpy.float32)

					# return
					return [nextDex,this_E,next_E,vlen,nextvlen,law,intt,this_T,this_U,this_Eedge,numpy.array([0]),numpy.array([0]),numpy.array([0])]

				else:  # return 0 if below the first energy]
					next_E = dataE[0]
					nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
					return [nextDex,0,0,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]

			elif law==4 or law == 44 or law==61:  
			# Kalbach-87 tabular distribution, or correlated angle-energy dist

				# get data
				if hasattr(rxn.energy_dist,"multiplicity"):
					m_yield    = rxn.energy_dist.multiplicity
				else:
					m_yield 	= 0
				dataE   = rxn.energy_dist.energy_in
				dataMu  = rxn.energy_dist.energy_out
				dataCDF = rxn.energy_dist.cdf
				dataPDF = rxn.energy_dist.pdf

				# check length
				assert dataE.__len__() > 0

				# find the index of the scattering table energy
				if this_E >= dataE[0] and this_E <= dataE[-1]:
					
					data_dex = numpy.where( dataE >= this_E )[0][0]

					#get energy of next bin
					if data_dex == dataE.__len__()-1:
						next_E  = self.MT_E_grid[-1]
						plusone = 0
					else:
						next_E  = dataE[data_dex+1]
						plusone = 1
					# find main E grid indext of next energy
					nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]

					# construct vector
					vlen  		=                         dataCDF[data_dex].__len__()
					cdf   		= numpy.ascontiguousarray(dataCDF[data_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
					pdf   		= numpy.ascontiguousarray(dataPDF[data_dex],dtype=numpy.float32)
					mu    		= numpy.ascontiguousarray(dataMu[ data_dex], dtype=numpy.float32)
					nextvlen 	=                         dataCDF[data_dex+ plusone].__len__()
					nextcdf  	= numpy.ascontiguousarray(dataCDF[data_dex+ plusone],dtype=numpy.float32) 
					nextpdf  	= numpy.ascontiguousarray(dataPDF[data_dex+ plusone],dtype=numpy.float32) 
					nextmu   	= numpy.ascontiguousarray(dataMu[ data_dex+ plusone],dtype=numpy.float32)
					intt 		= rxn.energy_dist.intt[data_dex]
					
					#check to make sure the same length
					assert vlen == mu.__len__()
					return [nextDex,this_E,next_E,vlen,nextvlen,law,intt,mu,cdf,pdf,nextmu,nextcdf,nextpdf]
				
				else:  # return 0 if below the first energy]
					next_E = dataE[0]
					nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
					return [nextDex,0,0,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
			elif law==66:
				#N-body phase space distribution

				# get data
				if hasattr(rxn.energy_dist,"multiplicity"):
					m_yield    = rxn.energy_dist.multiplicity
				else:
					m_yield 	= 0
				dataE   = rxn.energy_dist.energy
				intt    = 0

				# check length
				assert dataE.__len__() > 0

				# find the index of the scattering table energy
				if this_E >= dataE[0] and this_E <= dataE[-1]:
					
					#get energy of next bin		
					next_E   = self.MT_E_grid[-1]
					vlen     = 3.
					nextvlen = 3.

					# find main E grid indext of next energy
					nextDex = self.MT_E_grid.__len__()

					# data
					Q 			= numpy.ascontiguousarray(numpy.array([rxn.Q,rxn.energy_dist.nbodies,rxn.energy_dist.massratio]),dtype=numpy.float32)
										
					return [nextDex,this_E,next_E,vlen,nextvlen,law,intt,Q,Q,Q,Q,Q,Q]
				
				else:  # return 0 if below the first energy]
					next_E = dataE[0]
					nextDex = (numpy.where( self.MT_E_grid >= next_E )[0][0]) 
					return [nextDex,0,0,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
			else:
				print "law ",law," not handled, writing nulls"
				next_E = self.MT_E_grid[self.num_main_E-1]
				return [(self.MT_E_grid.__len__()-1),this_E,next_E,0,0,law,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
		else:
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has no energy tables"
			next_E   = self.MT_E_grid[self.num_main_E-1]
			nextDex = self.MT_E_grid.__len__()
			return [nextDex,this_E,next_E,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]



