#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void make_mask_kernel(unsigned N, unsigned * mask, unsigned * rxn,  unsigned begin, unsigned end){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	unsigned this_rxn = rxn[tid];

	if(this_rxn >= begin & this_rxn <= end){
		mask[tid] = 1;
	}
	else{
		mask[tid] = 0;
	}

}

void make_mask( unsigned NUM_THREADS, unsigned N, unsigned * mask, unsigned * rxn, unsigned begin, unsigned end){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	make_mask_kernel <<< blks, NUM_THREADS >>> (  N, mask, rxn, begin,  end);
	cudaThreadSynchronize();

}

