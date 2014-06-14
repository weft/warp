#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "LCRNG.cuh"

__global__ void rebase_yield_kernel(unsigned N, float keff, unsigned* rn_bank, unsigned* yield){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}
	if (yield[tid]==0){return;}
	//if(done[tid]){return;}

	unsigned 	this_yield 	=	yield[tid];
	unsigned 	rn 			=	rn_bank[tid];
	float 		new_yield 	=	(float)this_yield / keff;
	unsigned 	i_new_yield	=	(unsigned) new_yield;
	float 		rn1 		=	get_rand(&rn);

	if((float)i_new_yield+rn1 < new_yield){
		this_yield = i_new_yield+1;
	}
	else{
		this_yield = i_new_yield;
	}
	//printf("%u %6.4E %6.4E %u %u %6.4E\n",yield[tid],keff,new_yield,i_new_yield,this_yield,rn1);

	yield[tid] 	= 	this_yield;
	rn_bank[tid]= 	rn;

}

void rebase_yield( unsigned NUM_THREADS, unsigned N, float keff, unsigned* rn_bank, unsigned* yield){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	rebase_yield_kernel <<< blks, NUM_THREADS >>> (  N, keff, rn_bank, yield);
	cudaThreadSynchronize();

}

