#ifndef WARP_CUDA_H
#define WARP_CUDA_H

#include "device_copies.h"

// host calls
void print_banner();
void write_to_file(unsigned*,unsigned,std::string);

// device calls
void set_positions_rand( unsigned, unsigned, unsigned, source_point * , unsigned *  , float  * );
void copy_points( unsigned , unsigned , unsigned*  , unsigned  , unsigned *  , source_point *  , source_point * , float*, float*);
void sample_fission_spectra(unsigned , unsigned , unsigned , unsigned* , unsigned* , unsigned* , unsigned* , float * , float *, source_point* , float** );
void sample_fixed_source( unsigned,unsigned,unsigned*,unsigned*,float*,source_point*);
void macroscopic( unsigned , unsigned, unsigned, unsigned, unsigned, unsigned , unsigned*, source_point * , unsigned* , unsigned * , unsigned*, unsigned * , unsigned*, float * , unsigned * , float * , float *  , float* , unsigned*);
void microscopic( unsigned , unsigned, unsigned , unsigned , unsigned*, unsigned* , unsigned * , float * , unsigned * , float * , float *  , unsigned * , unsigned * ,  float* , unsigned * , float*, unsigned* );
void tally_spec( unsigned ,  unsigned, unsigned , unsigned , unsigned*, source_point * , float* , float* , float * , unsigned * , unsigned*, unsigned*, unsigned*, float*);
void escatter( cudaStream_t, unsigned , unsigned, unsigned , unsigned*, unsigned* , unsigned* , unsigned* , float*, source_point* , unsigned*, float*, float*, unsigned*, float**);
void iscatter( cudaStream_t, unsigned , unsigned , unsigned , unsigned*, unsigned* , unsigned * , unsigned * , float *, source_point *  ,unsigned * , float* , float* , unsigned* , float**, float**);
void cscatter( cudaStream_t, unsigned , unsigned, unsigned, unsigned , unsigned*, unsigned* , unsigned * , unsigned * , float *, source_point *  ,unsigned * , float* , float* , unsigned* , float**, float**);
void fission(  cudaStream_t , unsigned , unsigned , unsigned  , unsigned* , unsigned* , unsigned* , unsigned* , float* , source_point* , unsigned* , float* , unsigned* , float* , float** , float** );
void absorb(   cudaStream_t, unsigned , unsigned , unsigned*, unsigned*  , unsigned* );
void find_E_grid_index( unsigned , unsigned, unsigned , unsigned* , float * , float* , unsigned *, unsigned* );
void print_histories(unsigned, unsigned, unsigned *, unsigned*, source_point*, float*, unsigned*,unsigned*,unsigned*);
void pop_secondaries( unsigned, unsigned, unsigned* , unsigned* , unsigned* , unsigned* , unsigned*, unsigned*, source_point* , float* , unsigned* , float** );
void pop_source( unsigned ,  unsigned , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , source_point* , float*  , unsigned* , float ** , float** , source_point* , float* , float * , float* );
void rebase_yield( unsigned , unsigned , float , unsigned* , unsigned* );
void reaction_edges( unsigned ,  unsigned , unsigned* , unsigned* );
void check_remap( unsigned , unsigned , unsigned* , unsigned* , unsigned* );
void print_data( cudaStream_t , unsigned , unsigned , source_point* , float* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* , unsigned* );


#endif
