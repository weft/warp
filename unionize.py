from pyne import ace
from pyne import nucname
import numpy
import sys
import glob
import pylab

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
		## library temperature extension
		self.temp_extension   = '.03c'
		## cross section tables
		self.tables           = []
		## cross section libraries
		self.libraries        = []
		## AWR array
		self.awr 	      = []
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
	def _read_tables(self):

		datapath = '/usr/local/SERPENT/xsdata/endfb7/acedata/'
		
		for tope in self.isotope_list:
			#tope_number = nucname.mcnp(tope)
			#print tope
			#print nucname.mcnp(tope)
			#print glob.glob(datapath+str(tope_number)+'[A-Z]*[0-9]*.ace')
			librarypath=glob.glob(datapath+str(tope)+'[A-Z]*[0-9]*.ace')[0]
			self.libraries.append(ace.Library(librarypath))

		for lib in self.libraries:
			lib.read()
			iname=lib.tables.keys()[0][0:-4]   #strip off temp to get isotope name
			print "  loading "+iname+self.temp_extension
			self.tables.append(lib.find_table(iname+self.temp_extension))

		self.num_isotopes=self.libraries.__len__()

	##
	# \brief unionization function
	# \details unionizes MT energy grid and scattering energies in if present.
	# @param[in] self - material with attributes to be unionized
	def _unionize(self):

		for table in self.tables:
			self.MT_E_grid=numpy.union1d(self.MT_E_grid,table.energy)
			# unionize the scattering energies in as well!  if present of course
			for MT in table.reactions:
				rxn = table.reactions[MT]
				if hasattr(rxn,"ang_energy_in"):
					self.MT_E_grid=numpy.union1d(self.MT_E_grid,rxn.ang_energy_in)
				if hasattr(rxn,"energy_dist") and rxn.energy_dist.law!=3:
					self.MT_E_grid=numpy.union1d(self.MT_E_grid,rxn.energy_dist.energy_in)

		self.num_main_E   = self.MT_E_grid.__len__()
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

			#print "interpolating isotope "+str(tope_index)

			#do this isotopes entry in the total block
			this_array = numpy.interp( self.MT_E_grid, table.energy, table.sigma_t , left=0.0 )
			self.MT_array[:,tope_index]=this_array

			for MT in table.reactions:
				rxn        = table.reactions[MT]
				IE         = rxn.IE-1       #convert to python/C indexing 
				#print table.energy[IE:]
				#print rxn.sigma
				#if hasattr(rxn,'ang_energy_in'): 
				#	print rxn.ang_energy_in
				#else:
				#	print "no angular"
				#print rxn.threshold()
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
		#print "isotope = "+ str(isotope)+" = "+self.tables[isotope].name+" MT = "+str(MTnum)+" row="+str(row)+" col="+str(col)
		#print "energy="+str(self.MT_E_grid[row])
		if hasattr(rxn,"ang_energy_in"):
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has scattering data"
			scatterE   = rxn.ang_energy_in
			scatterMu  = rxn.ang_cos 
			scatterCDF = rxn.ang_cdf 
			# check length
			assert scatterE.__len__() > 0
			# find the index of the scattering table energy
			if this_E >= scatterE[0] and this_E <= scatterE[-1]:
				#print numpy.where( scatterE >= this_E )
				#print scatterE
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
				#print row,col
				#print this_E,next_E
				#print scatter_dex
				#print nextDex
				#print scatterMu [scatter_dex]
				#print scatterCDF[scatter_dex]
				# construct vector
				vlen      = scatterCDF[scatter_dex].__len__()
				cdf       = numpy.ascontiguousarray(scatterCDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
				mu        = numpy.ascontiguousarray(scatterMu[scatter_dex], dtype=numpy.float32)
				nextvlen  = scatterCDF[scatter_dex+plusone].__len__()
				nextcdf   = numpy.ascontiguousarray(scatterCDF[scatter_dex+ plusone],dtype=numpy.float32) 
				nextmu    = numpy.ascontiguousarray(scatterMu[scatter_dex + plusone], dtype=numpy.float32)
				#check to make sure the same lengths
				assert vlen == mu.__len__()
				# return
				#print "vlen="+str(vlen)
				return [nextDex,this_E,next_E,vlen,nextvlen,mu,cdf,nextmu,nextcdf]
			else:  # return 0 if below the first energy]
				next_E = scatterE[0]
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
				#print "energy starts at dex "+str(nextDex)+", energy="+str(next_E)+","+str(self.MT_E_grid[nextDex])
				return [nextDex,this_E,next_E,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
		elif hasattr(rxn,"energy_dist") and hasattr(rxn.energy_dist,"ang"):
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has 44 data"
			scatterE   = rxn.energy_dist.energy_in
			scatterMu  = rxn.energy_dist.ang
			scatterCDF = rxn.energy_dist.frac 
			# check length
			assert scatterE.__len__() > 0
			# find the index of the scattering table energy
			if this_E >= scatterE[0] and this_E <= scatterE[-1]:
				#print numpy.where( scatterE >= this_E )
				#print scatterE
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
				#print row,col
				#print this_E,next_E
				#print scatter_dex
				#print nextDex
				#print scatterMu [scatter_dex]
				#print scatterCDF[scatter_dex]
				# construct vector
				vlen      = scatterCDF[scatter_dex].__len__()
				cdf       = numpy.ascontiguousarray(scatterCDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
				mu        = numpy.ascontiguousarray(scatterMu[scatter_dex], dtype=numpy.float32)
				nextvlen  = scatterCDF[scatter_dex+plusone].__len__()
				nextcdf   = numpy.ascontiguousarray(scatterCDF[scatter_dex+ plusone],dtype=numpy.float32) 
				nextmu    = numpy.ascontiguousarray(scatterMu[scatter_dex + plusone], dtype=numpy.float32)
				#check to make sure the same lengths
				assert vlen == mu.__len__()
				# return
				#print "vlen="+str(vlen)
				#if MTnum == 91 :
				#	print "MT=91  thisE="+str(this_E)+" nextE="+str(next_E)+" thisdex="+str(row)+" nextdex="+str(nextDex)
				#	pylab.plot(mu,cdf)
				#	pylab.savefig("MT91_scatter_"+str(this_E)+".png")
				#	pylab.cla()
				return [nextDex,this_E,next_E,vlen,nextvlen,mu,cdf,nextmu,nextcdf]
			else:  # return 0 if below the first energy]
				next_E = scatterE[0]
				nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
				#if MTnum==91:
					#print "energy starts at dex "+str(nextDex)+", energy="+str(next_E)+","+str(self.MT_E_grid[nextDex])
				return [nextDex,this_E,next_E,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
		elif hasattr(table,"nu_t_energy"):
			# return interpolated nu values
			interped_nu = numpy.interp( self.MT_E_grid, table.nu_t_energy, table.nu_t_value )   #
			interped_nu = numpy.ascontiguousarray(interped_nu, dtype=numpy.float32)
			#print interped_nu
			#print "nu for MT="+str(MTnum)
			return [-1,-1,-1,-1,-1,interped_nu,interped_nu,interped_nu,interped_nu]
		else:
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has no angular tables"
			next_E   = self.MT_E_grid[self.num_main_E-1]
			nextDex = self.MT_E_grid.__len__()
			return [nextDex,this_E,next_E,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]



	##
	# \brief gets table of energy data
	# \details table returned in form of [nextDex, length, mu, cdf]
	# @param[in] self - isotope
	# @param[in] row - point in energy grid
	# @param[in] col - MT number
	def _get_energy_data(self,row,col):
		# scatter table returned in this form
		# returns [nextDex, length, mu, cdf] if scattering data exists

		#find the isotope we are in
		numbers = numpy.cumsum(self.reaction_numbers_total)
		isotope = 0
		for n in numbers:
			if (col - self.num_isotopes) <= n:
				break
			else:
				isotope = isotope + 1

		table = self.tables[isotope]
		MTnum = self.reaction_numbers[col]
		rxn   = table.reactions[MTnum]
		# get the energy from this index
		this_E = self.MT_E_grid[row]
		if hasattr(rxn,"energy_dist"):
			#print "LAW="+str(rxn.energy_dist.law)+" MT="+str(MTnum)
			if rxn.energy_dist.law == 3:
				next_E = self.MT_E_grid[self.num_main_E-1]
				return [(self.MT_E_grid.__len__()-1),this_E,next_E,0,0,3,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
			else:
				scatterE   = rxn.energy_dist.energy_in
				scatterMu  = rxn.energy_dist.energy_out
				scatterCDF = rxn.energy_dist.cdf
				scatterPDF = rxn.energy_dist.pdf
				law        = rxn.energy_dist.law
				#print "MT "+str(MTnum)+" law "+str(law)
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
					#print "isotope = "+ str(isotope)+" MT = "+str(MTnum)
					#print row,col
					#print this_E,next_E
					#print scatter_dex
					#print nextDex
					#print scatterMu [scatter_dex]
					#print scatterCDF[scatter_dex]
					# construct vector
					vlen  		= scatterCDF[scatter_dex].__len__()
					cdf   		= numpy.ascontiguousarray(scatterCDF[scatter_dex],dtype=numpy.float32)  # C/F order doesn't matter for 1d arrays
					pdf   		= numpy.ascontiguousarray(scatterPDF[scatter_dex],dtype=numpy.float32)
					mu    		= numpy.ascontiguousarray(scatterMu[scatter_dex], dtype=numpy.float32)
					nextvlen 	= scatterCDF[scatter_dex+plusone].__len__()
					nextcdf  	= numpy.ascontiguousarray(scatterCDF[scatter_dex+ plusone],dtype=numpy.float32) 
					nextpdf  	= numpy.ascontiguousarray(scatterPDF[scatter_dex+ plusone],dtype=numpy.float32) 
					nextmu   	= numpy.ascontiguousarray(scatterMu[scatter_dex + plusone], dtype=numpy.float32)
					#check to make sure the same length
					assert vlen == mu.__len__()
					#print "vlen,next "+str(vlen)+" "+str(nextvlen)
					# return
					#if MTnum == 91 :
					#	print "MT=91  thisE="+str(this_E)+" nextE="+str(next_E)+" thisdex="+str(row)+" nextdex="+str(nextDex)
					#	print scatter_dex
					#	pylab.plot(mu,cdf)
					#	pylab.savefig("MT91_energy_"+str(this_E)+".png")
					#	pylab.cla()
					return [nextDex,this_E,next_E,vlen,nextvlen,law,mu,cdf,pdf,nextmu,nextcdf,nextpdf]
				else:  # return 0 if below the first energy]
					next_E = scatterE[0]
					nextDex = numpy.where( self.MT_E_grid == next_E )[0][0]
					#if MTnum==91:
						#print "energy starts at dex "+str(nextDex)+", energy="+str(next_E)+","+str(self.MT_E_grid[nextDex])
					return [nextDex,0,0,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]
		else:
			#print "isotope "+str(isotope)+", MT = "+str(MTnum)+" has no energy tables"
			next_E   = self.MT_E_grid[self.num_main_E-1]
			nextDex = self.MT_E_grid.__len__()
			return [nextDex,this_E,next_E,0,0,0,numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0]),numpy.array([0])]



