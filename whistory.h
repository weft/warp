#ifndef WHISTORY_H
#define WHISTORY_H

/**
 * \class whistory whistory.h
 * \brief whistory class
 */

class whistory { 
	
	// wgeometry instance
	wgeometry					problem_geom;			/**< problem geometry */
	std::string					accel_type;				/**< acceleration type */
	
	// CUDPP
	CUDPPHandle					theCudpp;				/**< CUDPP handle */
	CUDPPConfiguration			compact_config;			/**< CUDPP compact configuration */
	CUDPPConfiguration			scan_int_config;		/**< CUDPP scan int configuration */
	CUDPPConfiguration			redu_int_config;		/**< CUDPP reduced int configuration */
	CUDPPConfiguration			redu_float_config; 		/**< CUDPP reduced float configuration */
	CUDPPConfiguration			radix_config;			/**< CUDPP radix configuration */
	CUDPPHandle					scanplan_int;			/**< CUDPP scan plan int handle */
	CUDPPHandle					reduplan_int;			/**< CUDPP reduce plan int handle */
	CUDPPHandle					reduplan_float;			/**< CUDPP reduce plan float handle */
	CUDPPHandle					compactplan;			/**< CUDPP compact plan handle */
	CUDPPHandle					radixplan;				/**< CUDPP radix plan handle */
	CUDPPResult					res;					/**< CUDPP result */
	
	// CURAND generator handle
	curandGenerator_t			rand_gen;				/**< random number generator handle */
	
	// cuda parameters
	unsigned					N;						/**< number of histories */
	unsigned					Ndataset;				/**< dataset size for number of histories */
	unsigned					NUM_THREADS;			/**< number of threads per block */
	unsigned					blks;					/**< number of blocks */
	unsigned					compute_device;			/**< compute device */
	cudaStream_t				stream[5];				/**< CUDA streams cor concurrent kernels */
	
	// host/device copied data
	cross_section_data*			d_xsdata;				/**< device cross section data structure containing device pointers */
	cross_section_data			dh_xsdata;				/**< host cross section data structure containing device pointers*/
	cross_section_data			h_xsdata;				/**< host cross section data structure containing host pointers*/
	particle_data*				d_particles;			/**< device particle data structure containing device pointers */
	particle_data				dh_particles;			/**< host particle data structure containing device pointers*/
	particle_data				h_particles;			/**< host particle data structure containing host pointers*/
	tally_data*					d_tally;				/**< device tally data structure containing device pointers*/
	tally_data*					dh_tally;				/**< host tally data structure containing device pointers*/
	tally_data*					h_tally;				/**< host tally data structure containing hist pointers*/

	// mapped arrays
	unsigned					n_edges;				/**< mapped array of number of edges */
	unsigned*					  edges;				/**< mapped array of edges */
	unsigned*					d_edges;				/**< device mapped array of edges */
	unsigned*					 reduced_yields;		/**< reduced yields */
	float*						 reduced_weight;		/**< reduced weight */
	unsigned*					d_reduced_yields;		/**< device reduced yields */
	float*						d_reduced_weight;		/**< device reduced weight */
	
	// 
    long unsigned				reduced_yields_total;   /**< long unsigned for accumulating yield numbers accurately on the host */
    double						reduced_weight_total;   /**< double for accumulating weight numbers accurately on the host */

	// materials data
	unsigned*					material_list;			/**< material list */
	unsigned*					isotope_list;			/**< isotope list */
	float*						number_density_matrix;	/**< isotope number density matrix */
	unsigned*					d_material_list;		/**< device material list */
	unsigned*					d_isotope_list;			/**< device isotope list */
	float*						d_number_density_matrix;/**< device isotope number density matrix */

	// reference remapping arrays
	unsigned*					  remap;				/**< remap */
	unsigned*					d_remap;				/**< remap pointer */

	// [re]initialization arrays
	unsigned*					zeros;					/**< zeros array */
	unsigned*					d_zeros;				/**< device zeros array */
	unsigned*					ones;					/**< int ones array */
	float*						fones;					/**< float ones array */

	// host-only data
	PyObject* 					xsdat_instance;			/**< Python object that loads and manipulates the cross section data */
	unsigned					RUN_FLAG;				/**< run flag */
	unsigned					outer_cell;				/**< outermost cell */
	unsigned					outer_cell_type;		/**< outermost cell type*/
	unsigned					n_materials;			/**< number of materials */
	unsigned					n_isotopes;				/**< number of isotopes */
	unsigned					n_tallies;				/**< number of tallies */
	unsigned					n_skip;					/**< number of cycles to skip */
	unsigned					n_cycles;				/**< number of active cycles */
	float						keff_sum;				/**< keff sum */
	float						keff2_sum;				/**< keff squared sum */
	float						keff_err;				/**< keff error */
	std::string					filename;				/**< file name */
	unsigned					is_initialized;			/**< init flag */
	unsigned					print_flag;				/**< print verbosity level*/
	unsigned					dump_flag;				/**< dump level*/

	// device-only variables
	unsigned*					d_valid_result;			/**< valid result pointer */
	unsigned*					d_valid_N;				/**< valied number of histories pointer */
	float*						d_fissile_energy;		/**< device fissile energy */
	spatial_data*				d_fissile_points;		/**< device fissile points */
 	unsigned*					d_scanned;				/**< device scanned pointer */
 	unsigned*					d_num_completed;		/**< device number of completed histories */
 	unsigned*					d_num_active;			/**< device number of active histories */
 	spatial_data*				d_bank_space;			/**< device bank space */
 	float*						d_bank_E;				/**< device bank energy */

	// xs data parameters used in parsing, copying, etc
	std::vector<std::string>	isotopes;				/**< cross section isotope string */
	std::vector<unsigned>		xs_num_rxns;			/**< cross section number of reactions */
	std::vector<unsigned>		xs_isotope_ints;		/**< cross section isotope numbers */
	unsigned					total_bytes_scatter;	/**< total size of scattering data */
	unsigned					total_bytes_energy ;	/**< total size of energy data */
	unsigned					MT_rows;				/**< MT number rows */
	unsigned					MT_columns;				/**< MT number columns */

	//geom parameters
	float						outer_cell_dims [6];	/**< outer cell minima and maxima */
	long unsigned*				fiss_img;				/**< fissile image accumulation */

	// private transport functions 
	/**
	 * \brief initializes the random number generator
	 */
	void init_RNG();
	/**
	 * \brief updates the random number
	 */
	void update_RNG();
	/**
	 * \brief initializes CUDPP
	 * \details initializes global objects, sorting stuff, int reduction stuff, 
	 * float reduction stuff, int scan stuff, radix sort stuff.
	 */
	void init_CUDPP();
	/**
	 * \brief initializes data on the host
	 * \details prepares data arrays by filling them with zeros
	 */
	void init_host();
		/**
	 * \brief initializes data on the device
	 * \details prepares data arrays by filling them with zeros
	 */
	void init_device();
	/**
	 * \brief copies data from the host device to the compute device
	 * \details copies history data, cross section data, and the device pointer
	 * array. zeros out the tally arrays.
	 */
	void copy_data_to_device();
	/**
	 * \brief loads cross sections
	 * \details makes isotope list, initializes the cross section libraries, reads
	 * the cross section tables, unionizes the main energy grid across all isotopes,
	 * makes the total MT reaction list from all isotopes, allocates the unionized
	 * array, inserts and interpolates the cross sections, gets the MT array buffer,
	 * gets and copies the unionized MT array, gets the unionized main energy grid
	 * buffer, gets the MT number vector, gets the MT number total vector, gets 
	 * the lengths vector, gets the AWR vector, gets the Q vector. does scattering
	 * stuff and energy stuff. passes information to the geometry.
	 */
	void init_cross_sections();
	/**
	 * \brief does an OptiX trace
	 * @param[in] type - trace type
	 */
	void trace(unsigned);
	/**
	 * \brief does an OptiX trace for a given number of active histories
	 * @param[in] type - trace type
	 * @param[in] n_active - number of active histories
	 */
	void trace(unsigned,unsigned);
	/**
	 * \brief reduces done values
	 * \returns reduced_done - number of done values
	 */
	unsigned reduce_done();
	/**
	 * \brief resets the cycle in criticality mode
	 * \details rebases the yield so that keff is 1, scans the yield to see where
	 * to write, sorts the reaction vector, populates the source, resets the run
	 * arrays, and updates the random numbers.
	 * @param[in] keff_cycle - keff value of previous cycle
	 */
	void reset_cycle(float);
	/**
	 * \brief resets the cycle in fixed-source mode
	 * \details resets the read-in run arrays, samples the fixed source, and 
	 * updates the random numbers.
	 */
	void reset_fixed();
	/**
	 * \brief not called in whistory.cpp
	 */
	void converge(unsigned);
	/**
	 * \brief samples fissile points
	 * \details updates the random numbers, sets uniformly random positions, runs
	 * OptiX to get the cell number, compacts data, copies the data back, copies
	 * new values in, writes starting positions to file, and updates the random
	 * numbers.
	 */
	void sample_fissile_points();
	/**
	 * \brief reduces yield values
	 * \returns total
	 */
	unsigned reduce_yield();
    /**
	 * \brief reduces weight values
	 * \returns total
	 */
    float reduce_weight();
    /**
     * \brief accumulates yields into host side values
     * @param[in] iteration - the active iteration number (starts at 0)
     * @param[in] keff - the running cumulative keff
     * @param[in] keff_cycle - the keff of the last cycle, used to renormalize the source vector
     */
    void accumulate_keff(unsigned, unsigned, double*, float*);
    /**
     * \brief raccumulates the flux tally
     * @param[in] 
     */
    void accumulate_tally();
	/**
	 * \brief returns how long it takes to do something
	 */
	float get_time();
	/**
	 * \brief prepares for secondary neutrons
	 * \details scans yields to determine where threads write into the done data,
	 * compacts done data to know where to write. 
	 */
	void prep_secondaries();
	/**
	 * \brief maps done histories
	 * \details flips done flag, remaps data to active histories, flips done flag 
	 * back.
	 * \returns num_active - number of active histories
	 */
	unsigned map_active();
	/**
	 * \brief remaps active histories
	 * \details sorts reaction values, launches edge detection kernels, calculates
	 * values for reaction blocks, calculates the total number of active histories,
	 * and rezeros the edge vector.
	 * @param[in] num_active - number of active histories
	 * @param[in] escatter_N - number of elastic scatters
	 * @param[in] escatter_start - elastic scattering start
	 * @param[in] iscatter_N - number of inelastic scatters
	 * @param[in] iscatter_start - inelastic scattering start
	 * @param[in] cscatter_N - number of compound scatters
	 * @param[in] cscatter_start - compound scattering start
	 * @param[in] fission_N - number of fissions
	 * @param[in] fission_start - fission start
	 */
	void  remap_active(unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* );
	/**
	 * \brief prints the locations of the source points to file
	 * @param[in] array_in - source point array
	 * @param[in] N - number of histories
	 * @param[in] filename - filename
	 * @param[in] opentype - file extension
	 */
	void  write_to_file(spatial_data*  , unsigned , std::string , std::string);
	/**
	 * \brief prints the locations of the source points to file
	 * @param[in] array_in - source point array
	 * @param[in] array_in2 - second array
	 * @param[in] N - number of histories
	 * @param[in] filename - filename
	 * @param[in] opentype - file extension
	 */
	void  write_to_file(spatial_data*  , float*, unsigned , std::string , std::string);
	/**
	 * \brief prints the source points to a file
	 * @param[in] array_in - source point array
	 * @param[in] N - number of histories
	 * @param[in] filename - filename
	 * @param[in] opentype - file extension
	 */
	void  write_to_file(unsigned*  , unsigned , std::string, std::string );
	/**
	 * \brief prints the locations of the source points to file
	 * @param[in] array_in - source point array
	 * @param[in] array_in2 - second array
	 * @param[in] N - number of histories
	 * @param[in] filename - filename
	 * @param[in] opentype - file extension
	 */
	void  write_to_file(unsigned*  , unsigned*, unsigned , std::string, std::string );
	/**
	 * \brief writes results to file
	 * @param[in] runtime - runtime
	 * @param[in] keff - keff
	 * @param[in] opentype - file extension
	 */
	void  write_results(float,float,std::string);
	/**
	 * \brief calls python function, copys returned buffer to C and CUDA pointers
	 * @param[in] 
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_python_buffer(float**,float**,std::string);
	/**
	 * \brief calls python function, copys returned buffer to C and CUDA pointers
	 * @param[in] 
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_python_buffer(unsigned**,unsigned**,std::string);
	/**
	 * \brief calls python function, copys returned buffer to C pointer (no cuda)
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_python_buffer(float**,std::string);
	/**
	 * \brief calls python function, copys returned buffer to C pointer (no cuda)
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_python_buffer(unsigned**,std::string);
	 /**
	 * \brief initialized cross section data object in python
	 * @param[in] 
	 * @param[out]
	 */
	int init_python();
	/**
	 * \brief calls python function, copys returned buffer to C and CUDA pointers
	 * @param[in] 
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_scatter_data();
	/**
	 * \brief calls python function, copys returned buffer to C and CUDA pointers
	 * @param[in] 
	 * @param[in] 
	 * @param[in] 
	 */
	void copy_energy_data();
public:
	/**
	* \brief constructor
	* \details makes geometry, sets tally vector length, creates dataset size, sets
	* compute device and acceleration type, creates CUDA streams.
	*/
	whistory(unsigned,wgeometry);
	/**
	* \brief destructor
	*/
	~whistory();
	/**
	* \brief prints cross section data information
	*/
	void print_xs_data();
	/**
	* \brief prints pointer information
	*/
	void print_pointers();
	/**
	 * \brief prints table of properties of geometry materials
	 */
	void print_materials_table();
	/**
	* \brief runs history
	* \details initializes run variables, clears fissile points file, opens run
	* stats file, records stats. finds the material and nearest surfact distance,
	* finds the main energy grid index, finds interaction length, computes spectra,
	* finds reaction type, remaps threads, does scattering reactions, reduces the
	* yield, resets the cycle, recalculates the running average, and prints the
	* total transport runtime.
	*/
	void run();
	/**
	* \brief writes cross section data to file
	* @param[in] filename - filename
	*/
	void write_xs_data(std::string);
	/**
	* \brief writes tally values to file
	* @param[in] tallynum - tally number
	*/
	void write_tally(unsigned);
	/**
	* \brief sets tally cell to input value
	* @param[in] cell - tally cell
	*/
	void set_tally_cell(unsigned);
	/**
	* \brief sets run type
	* @param[in] type_in - run type
	*/
	void set_run_type(unsigned);
	/**
	* \brief sets run type
	* @param[in] type_in - run type
	*/
	void set_run_type(std::string);
	/**
	* \brief sets number of cycles to skip and number of active cycles
	* @param[in] n_cycles_in - number of active cycles
	* @param[in] n_skip_in - number of cycles to skip
	*/
	void set_run_param(unsigned,unsigned);
	/**
	* \brief initialization function
	* \details initializes OptiX stuff and CUDA stuff, allocates device data, 
	* creates host data arrays, initializes counters to zero, copies outermost
	* cell dimensions and isotope list, maps edge array, initializes host values,
	* initializes the random number generator and CUDPP, loads cross sections, and
	* copies data to the compute device.
	*/
	void init();
	/**
	* \brief prints out details (model, memory, compute capability, etc.) of all
	* available compute devices
	*/
	void device_report();
	/**
	* \brief sets device number to input value
	* @param[in] dev_in - device number
	*/
	void set_device(unsigned);
	/**
	* \brief does nothing
	* @param[in] accel_in - acceleration type
	*/
	void set_acceration(std::string);
	/**
	* \brief sets filename to input string
	* @param[in] filename_in - filename
	*/
	void set_filename(std::string);
	/**
	* \brief appends history data to file in debug mode.
	* @param[in] iteration - iteration number
	*/
	void write_histories(unsigned iteration);
	/**
	* \brief sets amount of information printed to stdout
	* @param[in] level - level of verbosity
	*/
	void set_print_level(unsigned level);
	/**
	* \brief sets what types of information are dumped to files
	* @param[in] level - dump type flag
	*/
	void set_dump_level(unsigned level);
	/**
	* \brief produces png images of the geometry, named filename-[xy,xz,yz].png
	* @param[in] type, color is based 'cell' or 'material'
	*/
	void plot_geom(std::string type);
	/**
	* \brief creates a color map
	* @param[in] color - rgb colors, float[3]
	* @param[in] x - value 
	* @param[in] min,max - values used to normalize the color  
	*/
	void make_color(float* , unsigned , unsigned , unsigned );
	/**
	* \brief creates a hot2 color map
	* @param[in] color - rgb colors, float[3]
	* @param[in] x - value 
	* @param[in] min,max - values used to normalize the color  
	*/
	void hot2(float* , long unsigned , long unsigned , long unsigned );
	/**
	* \brief creates a binary colormap, black iff 0
	* @param[in] color - rgb colors, float[3]
	* @param[in] x - value 
	* @param[in] min,max - values used to normalize the color; unused, only present to keep arguments the same as other colormaps 
	*/
	void nonzero(float* , unsigned , unsigned , unsigned );
	/**
	* \brief bins and accumulates fission points to grid
	* @param[in] d_space - device space points
	* @param[in] N - dataset size 
	*/
	void bin_fission_points( spatial_data * , unsigned );
	/**
	* \brief writes binned fission point image to a .png
	*/
	void write_fission_points();
	/**
	* \brief prints an amazing and beautiful WARP banner to stdout
	*/
	void print_banner();
	/**
	* \brief does size logic on python buffer and writes values into pointers passed.  Made since python interface routines don't set first/second shape value if the array size is 0/1
	*/
	void get_Py_buffer_dims(unsigned* , unsigned* , unsigned* , Py_buffer* );
};

#endif
