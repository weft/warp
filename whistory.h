#ifndef WHISTORY_H
#define WHISTORY_H
/**
 * \class whistory whistory.h
 * \brief whistory class
 */

class whistory { 
	wgeometry 	       problem_geom; /**< problem geometry */
	std::string 	       accel_type;   /**< acceleration type */
	// CUDPP
	CUDPPHandle            theCudpp; /**< CUDPP handle */
	CUDPPHashTableConfig   hash_config; /**< CUDPP hash table configuration */
	CUDPPConfiguration     compact_config;    /**< CUDPP compact configuration */
	CUDPPConfiguration     scan_int_config;   /**< CUDPP scan int configuration */
	CUDPPConfiguration     redu_int_config;   /**< CUDPP reduced int configuration */
	CUDPPConfiguration     redu_float_config; /**< CUDPP reduced float configuration */
	CUDPPConfiguration     radix_config;	       /**< CUDPP radix configuration */
	CUDPPHandle            mate_hash_table_handle; /**< CUDPP material hash table handle */
	CUDPPHandle            fiss_hash_table_handle; /**< CUDPP fissile hash table handle */
	CUDPPHandle            scanplan_int; /**< CUDPP scan plan int handle */
	CUDPPHandle            reduplan_int; /**< CUDPP reduce plan int handle */
	CUDPPHandle            reduplan_float; /**< CUDPP reduce plan float handle */
	CUDPPHandle            compactplan; /**< CUDPP compact plan handle */
	CUDPPHandle            radixplan; /**< CUDPP radix plan handle */
	CUDPPResult            res; /**< CUDPP result */
	unsigned * 	       d_valid_result; /**< valid result pointer */
	unsigned * 	       d_valid_N; /**< valied number of histories pointer */
	unsigned * 	       d_remap; /**< remap pointer */
	// CURAND
	curandGenerator_t rand_gen; /**< random number generator */
	// cuda parameters
	unsigned 	N; /**< number of histories */
	unsigned 	Ndataset; /**< dataset size for number of histories */
	unsigned  	RNUM_PER_THREAD; /**< random numbers per thread */
	unsigned 	NUM_THREADS; /**< number of threads */
	unsigned 	blks; /**< number of blocks */
	unsigned 	compute_device; /**< compute device (always 0) */
	cudaStream_t 	stream[5]; /**< CUDA stream */
	// host data
	unsigned 	RUN_FLAG; /**< run flag */
	unsigned 	qnodes_depth; /**< quaternary node depth */
	unsigned	qnodes_width; /**< quaternary node width */
	unsigned    outer_cell; /**< outermost cell */
	unsigned    outer_cell_type; /**< outermost cell type*/
	unsigned 	n_materials; /**< number of materials */
	unsigned 	n_isotopes; /**< number of isotopes */
	unsigned 	n_tally; /**< number of tallies */
	unsigned 	n_qnodes; /**< number of quaternary nodes */
	unsigned 	n_skip; /**< number of cycles to skip */
	unsigned 	n_cycles; /**< number of active cycles */
	float 		keff_sum;
	float 		keff2_sum;
	float 		keff_err;
	std::string     filename; /**< file name */
	unsigned 	is_initialized;  /**< init flag */
	source_point *  space; /**< source point spatial pointer */
	unsigned  	print_flag; /**< print verbosity level*/
	unsigned 	dump_flag; /**< dump level*/
	/**
	 * \brief cross section length numbers
	 * \details 0 = isotopes, 1 = main menergy points, 2 = total number of reaction channels, 
	 * 3 = matrix energy points, 4 = angular cosine points, 5 = outgoing energy points
	 */
	unsigned *      xs_length_numbers;     // 0=isotopes, 1=main E points, 2=total numer of reaction channels, 3=matrix E points, 4=angular cosine points, 5=outgoing energy points 
	unsigned *      xs_MT_numbers_total; /**< total cross sextion MT numbers */
	unsigned *   	xs_MT_numbers; /**< cross section MT numbers */
	float *		xs_data_MT; /**< cross section MT data */
	float *		xs_data_main_E_grid; /**< cross section data main energy grid */
	float **	xs_data_scatter; /**< scattering cross section data */
	float **	xs_data_energy; /**< energy cross section data */
	float **	xs_data_scatter_host; /**< scattering cross section host data */
	float **	xs_data_energy_host; /**< energy cross section host data */
	float * 	xs_data_Q; /**< cross section data Q-values */
	float *         E; /**< energy */
	float *         Q; /**< Q-value */
	unsigned *	rn_bank; /**< reaction bank */ 
	float * 	awr_list; /**< atomic weight ratio (AWR) list */
	float *          tally_score;        /**< tally score */
	float *          tally_square;       /**< tally square */
	unsigned *       tally_count;        /**< tally count */
    double *         tally_score_total;  /**< tally score accumulated total */
    double *         tally_square_total; /**< tally square accumulaed total */
    long unsigned *  tally_count_total;  /**< tally count accumulated total */
	unsigned 	 tally_cell; /**<tally cell */
	unsigned * 	index; /**< index */
	unsigned *      cellnum; /**< cell number */
	unsigned *      matnum; /**< material number */
	unsigned *      isonum; /**< isotope number */
	unsigned *      rxn; /**< reaction */
	unsigned *      done; /**< done flag */
	unsigned *      yield; /**< yield */
	//unsigned * 	material_list; /**< material list */
	//unsigned * 	isotope_list; /**< isotope list */
	float *  	number_density_matrix; /**< isotope number density matrix */
	unsigned 	reduced_yields; /**< reduced yields */
	unsigned * 	remap; /**< remap */
	unsigned * 	zeros; /**< zeros */
	unsigned * 	ones; /**< ones */
    long unsigned   reduced_yields_total;   /**< long unsigned for accumulating yield numbers accurately on the host */
	qnode * 	qnodes; /**< quaternary nodes */
	// device data
	source_point *  d_space; /**< device spatial source point */
	unsigned *      d_xs_length_numbers; /**< device cross section length numbers */
	unsigned * 	d_xs_MT_numbers_total; /**< device cross section total MT numbers */
	unsigned * 	d_xs_MT_numbers; /**< device cross section MT numbers */
	float *		d_xs_data_MT; /**< device cross section MT numbers */
	float *		d_xs_data_main_E_grid; /**< device cross section main energy grid */
	float **	d_xs_data_scatter; /**< device scattering cross section data */
	float ** 	d_xs_data_energy; /**< device energy cross section data */
	float * 	d_xs_data_Q; /**< device cross section Q-value data */
	float *         d_E; /**< device energy */
	float *         d_Q; /**< device Q-value */
	unsigned *      d_rn_bank; /**< device reaction bank */
	float * 	d_awr_list; /**< device AWR list */
	float * 	d_tally_score; /**< device tally score */
	float *         d_tally_square; /**< device tally square */
	unsigned * 	d_tally_count; /**< device tally count */
	unsigned * 	d_index; /**< device index */
	unsigned *      d_cellnum; /**< device cell number */
	unsigned *      d_matnum; /**< device material number */
	unsigned *      d_isonum; /**< device isotope number */
	unsigned *      d_rxn; /**< device reaction */
	unsigned *      d_rxn_remap; /**< device reaction remap */
	unsigned *      d_done; /**< device done flag */
	unsigned *      d_yield; /**< device yield */
	unsigned * 	d_material_list; /**< device material list */
	unsigned * 	d_isotope_list; /**< device isotope list */
	float *  	d_number_density_matrix; /**< device isotope number density matrix */
	unsigned * 	d_reduced_yields; /**< device reduced yields */
	unsigned * 	d_reduced_done; /**< device reduced done flags */
	float * 	d_fissile_energy; /**< device fissile energy */
	source_point * 	d_fissile_points; /**< device fissile points */
	unsigned * 	d_mask; /**< device mask */
	qnode *		d_qnodes_root; /**< device quaternary nodes root */
	unsigned * 	d_completed; /**< device completed pointer */
 	unsigned *  	d_scanned; /**< device scanned pointer */
 	unsigned * 	d_active; /**< device active pointer */
 	unsigned * 	d_num_completed; /**< device number of completed histories */
 	unsigned * 	d_num_active; /**< device number of active histories */
 	source_point *  d_bank_space; /**< device bank space */
 	float * 	d_bank_E; /**< device bank energy */
    unsigned *  d_zeros; /**< zeros */
	// mapped arrays
	unsigned        n_edges; /**< mapped array of number of edges */
	unsigned*         edges; /**< mapped array of edges */
	unsigned*       d_edges; /**< device mapped array of edges */
	// xs data parameters
	std::string xs_isotope_string; /**< cross section isotope string */
	std::vector<unsigned> 	xs_num_rxns;     /**< cross section number of reactions */
	std::vector<unsigned> 	xs_isotope_ints; /**< cross section isotope numbers */
	unsigned 		total_bytes_scatter; /**< total size of scattering data */
	unsigned 		total_bytes_energy ; /**< total size of energy data */
	unsigned 		MT_rows;    /**< MT number rows */
	unsigned 		MT_columns; /**< MT number columns */
	//geom parameters
	float 		outer_cell_dims [6]; /**< outer cell minima and maxima */
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
	 * \brief initializes data on the host device
	 * \details prepares data arrays by filling them with zeros
	 */
	void init_host();
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
	void load_cross_sections();
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
	 * \brief builds a quad tree for energy search
	 */
	void create_quad_tree();
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
	void  write_to_file(source_point*  , unsigned , std::string , std::string);
	/**
	 * \brief prints the locations of the source points to file
	 * @param[in] array_in - source point array
	 * @param[in] array_in2 - second array
	 * @param[in] N - number of histories
	 * @param[in] filename - filename
	 * @param[in] opentype - file extension
	 */
	void  write_to_file(source_point*  , float*, unsigned , std::string , std::string);
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
};

#endif
