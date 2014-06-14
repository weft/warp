#include <cuda.h>
#include <stdio.h>

__global__ void flip_done_kernel(unsigned Ndataset, unsigned * done){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= Ndataset){return;}       //return if out of bounds

	unsigned bit = 1;

	if(done[tid]){bit=0;}

	done[tid]=bit;

}

void flip_done(unsigned NUM_THREADS, unsigned Ndataset, unsigned* d_done){

	unsigned blks = ( Ndataset + NUM_THREADS - 1 ) / NUM_THREADS;

	flip_done_kernel <<< blks, NUM_THREADS >>> ( Ndataset, d_done);
	cudaThreadSynchronize();

}

