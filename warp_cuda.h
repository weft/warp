#ifndef WARP_CUDA_H
#define WARP_CUDA_H

#include "device_copies.h"

////////////////
// host calls //
////////////////

/**
 * \brief prints warp banner
 * \details prints the WARP ASCII art banner to stdout... that's it.
 */
void print_banner();

/**
 * \brief writes CUDA array to a text file
 * \details copies the cuda array to a local buffer, writes to buffer to a new file, then frees the local memory
 * @param[in] array_in device pointer to array to write
 * @param[in] N        number of elements to write
 * @param[in] filename name for the file
 */ 
void write_to_file(unsigned*,unsigned,std::string);

///////////////////////////////
// device calls (c wrappers) //
///////////////////////////////

/**
 * \brief sets starting cycle points uniformly random
 * \details sets starting cycle points uniformly random in the specified volume
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    outer_cell_type the outer cell type, sets the shape of the sampling
 * @param[in]    d_space         device pointer to spatial data array
 * @param[in]    d_rn_bank       device pointer to random number array
 * @param[in]    outer_cell_dims host pointer to array of outer cell extrema
 */ 
void set_positions_rand( unsigned, unsigned, unsigned, spatial_data * , unsigned *  , float  * );

/**
 * \brief   copy points between two sets of space and energy data buffers, redirected with a mapping array
 * \details copy points between two sets of space and energy data buffers, redirected with a mapping array
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    Nout            the total number of threads to launch on the grid
 * @param[in]    Nvalid          the total number of device elements to copy from
 * @param[in]    current_index   starting index
 * @param[in]    to_valid        device pointer to data remapping vector
 * @param[in]    positions_out   device pointer to spatial data array destination
 * @param[in]    positions_in    device pointer to spatial data array source
 * @param[in]    E_out           device pointer to energy data array destination
 * @param[in]    E_in            device pointer to energy data array source
 */ 
void copy_points( unsigned , unsigned , unsigned*  , unsigned  , unsigned *  , spatial_data *  , spatial_data * , float*, float*);

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    active          device pointer to remapping vector
 * @param[in]    rn_bank         device pointer to random number array
 * @param[in]    E               device pointer to energy data array
 * @param[in]    space           device pointer to spatial data array
 */ 
void sample_fixed_source( unsigned,unsigned,unsigned*,unsigned*,float*,spatial_data*);

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS             the number of threads to run per thread block
 * @param[in]    N                       the total number of threads to launch on the grid
 * @param[in]    converged               flag for tally scoring
 * @param[in]    n_materials             number of materials
 * @param[in]    n_isotopes              number of isotopes
 * @param[in]    n_tallies               number of tallies
 * @param[in]    d_xsdata                device pointer to cross section data pointer array 
 * @param[in]    d_particles             device pointer to particle data pointer array 
 * @param[in]    d_tally                 device pointer to tally array
 * @param[in]    d_remap                 device pointer to data remapping vector
 * @param[in]    d_number_density_matrix device pointer to material number density array
 */ 
void macro_micro( unsigned , unsigned , unsigned ,  unsigned , unsigned , unsigned , cross_section_data* , particle_data* , tally_data* , unsigned* , float* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    stream          CUDA stream to launch the kernel on
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid for level scattering
 * @param[in]    starting_index  starting index of the level scatter block in the remap vector
 * @param[in]    d_xsdata        device pointer to cross section data pointer array 
 * @param[in]    d_particles     device pointer to particle data pointer array 
 * @param[in]    d_remap         device pointer to data remapping vector
 */ 
void scatter_level( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    stream          CUDA stream to launch the kernel on
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid for continuum scattering
 * @param[in]    starting_index  starting index of the continuum scatter block in the remap vector
 * @param[in]    d_xsdata        device pointer to cross section data pointer array 
 * @param[in]    d_particles     device pointer to particle data pointer array 
 * @param[in]    d_remap         device pointer to data remapping vector
 */ 
void scatter_conti( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    stream          CUDA stream to launch the kernel on
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid for multiplicity scattering
 * @param[in]    starting_index  starting index of the multiplicity scatter block in the remap vector
 * @param[in]    d_xsdata        device pointer to cross section data pointer array 
 * @param[in]    d_particles     device pointer to particle data pointer array 
 * @param[in]    d_remap         device pointer to data remapping vector
 */ 
void scatter_multi( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    stream          CUDA stream to launch the kernel on
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid for fission
 * @param[in]    starting_index  starting index of the fission block in the remap vector
 * @param[in]    d_xsdata        device pointer to cross section data pointer array 
 * @param[in]    d_particles     device pointer to particle data pointer array 
 * @param[in]    d_remap         device pointer to data remapping vector
 */ 
void fission(       cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    d_xsdata        device pointer to cross section data pointer array 
 * @param[in]    d_remap         device pointer to data remapping vector
 * @param[in]    d_E             device pointer to energy data array
 * @param[in]    d_index         device pointer to index array (stores the unionized grid index of the current energy)
 * @param[in]    d_rxn           device pointer of the reaction number array
 */ 
void find_E_grid_index( unsigned , unsigned, cross_section_data* , unsigned* , float* , unsigned *, unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS        the number of threads to run per thread block
 * @param[in]    N                  the total number of threads to launch on the grid
 * @param[in]    d_xsdata           device pointer to cross section data pointer array
 * @param[in]    d_particles        device pointer to particle data pointer array
 * @param[in]    d_scanned          device pointer to array of the cumulative sum (scan) of the yield array, used to find final index where new particles will be written
 * @param[in]    fission_particles  device pointer to intermadiate spatial data array where popped values will be written
 * @param[in]    fission_energy     device pointer to intermadiate energy data array where popped values will be written
 */ 
void pop_fission(  unsigned , unsigned , cross_section_data* , particle_data* , unsigned*, spatial_data*, float*) ;

/**
 * \brief 
 * \details 
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    keff            k-effective of the current cycle
 * @param[in]    d_particles     device pointer to particle data pointer array
 */ 
void rebase_yield( unsigned , unsigned , float , particle_data*);

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    d_edges         device pointer to the edges array (stores to indices of edges the reaction blocks)
 * @param[in]    d_rxn           device pointer of the reaction number array
 */ 
void reaction_edges( unsigned ,  unsigned , unsigned* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    a               Watt spectrum, parameter a
 * @param[in]    b               Watt spectrum, parameter b
 * @param[in]    rn_bank         device pointer to random number array
 * @param[in]    E               device pointer to energy data array
 */ 
void sample_fissile_energy(unsigned ,  unsigned , float , float , unsigned * , float*);

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    N               the total number of threads to launch on the grid
 * @param[in]    d_xsdata        device pointer to cross section data pointer array
 * @param[in]    d_particles     device pointer to particle data pointer array
 * @param[in]    d_remap         device pointer to data remapping vector
 */ 
void safety_check( unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS     the number of threads to run per thread block
 * @param[in]    dex0            starting index 
 * @param[in]    dex1            ending index
 * @param[in]    d_xsdata        device pointer to cross section data pointer array
 */ 
void check_pointers( unsigned, unsigned, unsigned, cross_section_data* );

#endif
