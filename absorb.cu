#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void absorb_kernel(unsigned N, unsigned* active, unsigned * rxn , unsigned* done){


	//PLACEHOLDER FOR FISSIONS, NEED TO READ NU TABLES LATER
	
	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}         //return if out of bounds
	
	//remap to active
	//tid=active[tid];
	if(done[tid]){return;}

	if (rxn[tid] < 102 | rxn[tid] > 800 ){return;}  //return if not some sort of absorption, ie (n,not-n)

	//printf("in abs, rxn=%u\n",rxn[tid]);

	done[tid]  = 1;

}

void absorb( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned* active, unsigned * rxn , unsigned* done){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	//absorb_kernel <<< blks, NUM_THREADS >>> (  N, active, rxn , done);
	absorb_kernel <<< blks, NUM_THREADS , 0 , stream >>> (  N, active, rxn , done);
	cudaThreadSynchronize();

}

