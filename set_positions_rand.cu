#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void set_positions_rand_kernel(unsigned N , unsigned RNUM_PER_THREAD, source_point * positions_ptr , unsigned * rn_bank , float x_min , float y_min , float z_min , float x_max , float y_max , float z_max ){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid>=N){return;}

	unsigned 	rn = rn_bank[tid];
	float rn1 = get_rand(&rn);
	float rn2 = get_rand(&rn);
	float rn3 = get_rand(&rn);
	float mu    = 2.0 * get_rand(&rn)-1.0;
	float theta = 6.28318530718 * get_rand(&rn);

	positions_ptr[tid].surf_dist =     500000;   
	positions_ptr[tid].x         =     0.9 * ( ( x_max - x_min ) * rn1 + x_min );  
	positions_ptr[tid].y         =     0.9 * ( ( y_max - y_min ) * rn2 + y_min );  
	positions_ptr[tid].z         =     0.9 * ( ( z_max - z_min ) * rn3 + z_min ); 
	positions_ptr[tid].xhat      =     sqrtf(1.0-mu*mu) * cosf( theta );
	positions_ptr[tid].yhat      =     sqrtf(1.0-mu*mu) * sinf( theta );
	positions_ptr[tid].zhat      =     mu;
  
	rn_bank[tid]	=	rn;

}

void set_positions_rand( unsigned NUM_THREADS, unsigned N, unsigned RNUM_PER_THREAD, source_point * d_space , unsigned * d_rn_bank, float * outer_cell_dims){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	set_positions_rand_kernel <<<  blks, NUM_THREADS >>> ( N , RNUM_PER_THREAD, d_space , d_rn_bank, outer_cell_dims[0], outer_cell_dims[1], outer_cell_dims[2], outer_cell_dims[3], outer_cell_dims[4], outer_cell_dims[5]);
	cudaThreadSynchronize();

}