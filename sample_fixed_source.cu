#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void sample_fixed_source_kernel(unsigned N, unsigned* active, unsigned* rn_bank , float * E, source_point* space){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	//remap to active
	//tid=active[tid];

	// load in
	unsigned rn   		=  rn_bank[ tid ];
	const float pi   	=   3.14159265359 ;
	const float mu   	= ( get_rand(&rn) ) * 2.0 - 1.0;
	const float theta	= ( get_rand(&rn) ) * 2.0 * pi ;

	//monoenergetic for now
	E[tid]=10.0;  //1.0e-6;

	//point source for now
	space[tid].x = 0.0;
	space[tid].y = 0.0;
	space[tid].z = 0.0;

	//set isotropic for now
	space[tid].xhat = sqrtf(1-mu*mu) * cosf( theta );
	space[tid].yhat = sqrtf(1-mu*mu) * sinf( theta );
    space[tid].zhat = mu;

    rn_bank[tid]	=	rn;


}

void sample_fixed_source( unsigned NUM_THREADS, unsigned N, unsigned* active, unsigned* rn_bank, float * E, source_point* space){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	sample_fixed_source_kernel <<< blks, NUM_THREADS >>> (  N, active, rn_bank, E , space );
	cudaThreadSynchronize();

}

