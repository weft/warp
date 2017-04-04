#ifndef WARP_CUDA_H
#define WARP_CUDA_H

#include "device_copies.h"

//
//
// ---> doxygen comments are done in the .cu files, doing it here never worked right
//
//

////////////////
// host calls //
////////////////


void write_to_file(unsigned*,unsigned,std::string);

///////////////////////////////
// device calls (c wrappers) //
///////////////////////////////


void set_positions_rand( unsigned, unsigned, unsigned, spatial_data * , unsigned *  , float  * );
void copy_points( unsigned , unsigned , unsigned*  , unsigned  , unsigned *  , spatial_data *  , spatial_data * , float*, float*);
void sample_fixed_source( unsigned,unsigned,unsigned*,unsigned*,float*,spatial_data*);
void macro_micro( unsigned , unsigned , unsigned ,  unsigned , unsigned , unsigned , cross_section_data* , particle_data* , tally_data* , unsigned* , float* );
void scatter_level( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );
void scatter_conti( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );
void scatter_multi( cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );
void fission(       cudaStream_t, unsigned , unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );
void find_E_grid_index( unsigned , unsigned, cross_section_data* , unsigned* , float* , unsigned *, unsigned* );
void pop_fission(  unsigned , unsigned , cross_section_data* , particle_data* , unsigned*, spatial_data*, float*) ;
void rebase_yield( unsigned , unsigned , float , particle_data*);
void reaction_edges( unsigned ,  unsigned , unsigned* , unsigned* );
void sample_fissile_energy(unsigned ,  unsigned , float , float , unsigned * , float*);
void safety_check( unsigned, unsigned, cross_section_data* , particle_data* , unsigned* );
void check_pointers( unsigned, unsigned, unsigned, cross_section_data* );

#endif
