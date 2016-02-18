#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"

__global__ void rebase_yield_kernel(unsigned N, float keff, particle_data* d_particles ){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// declare shared variables	
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			yield;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		rn_bank						= d_particles[0].rn_bank;
		yield						= d_particles[0].yield;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	if (yield[tid]==0){return;}

	unsigned 	this_yield 	=	yield[tid];
	unsigned 	rn 			=	rn_bank[tid];
	float 		new_yield 	=	(float) this_yield / keff;
	unsigned 	i_new_yield	=	(unsigned) new_yield;
	float 		rn1 		=	get_rand(&rn);

	if((float)i_new_yield+rn1 < new_yield){
		this_yield = i_new_yield+1;
	}
	else{
		this_yield = i_new_yield;
	}

	yield[tid] 	= 	this_yield;
	rn_bank[tid]= 	rn;

}

void rebase_yield( unsigned NUM_THREADS, unsigned N, float keff, particle_data* d_particles ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	rebase_yield_kernel <<< blks, NUM_THREADS >>> ( N, keff, d_particles );
	cudaThreadSynchronize();

}

