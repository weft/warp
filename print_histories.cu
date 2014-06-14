#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "LCRNG.cuh"

__global__ void print_histories_kernel( unsigned N , unsigned* isonum, unsigned* rxn, source_point* space, float* E, unsigned* done, unsigned* yield, unsigned* rn_bank){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}       //return if out of bounds
	//if (done[tid]){return;}
	unsigned rn = rn_bank[tid];
	
	printf("tid=%d isonum=%u rxn=%u yield=%u (x,y,z)= % 6.4E % 6.4E % 6.4E, E=% 6.4E, rn %u %6.4E\n",tid,isonum[tid],rxn[tid],yield[tid],space[tid].x,space[tid].y,space[tid].z,E[tid],rn,get_rand(&rn));

}

void print_histories( unsigned NUM_THREADS, unsigned N, unsigned* isonum, unsigned* rxn, source_point* space, float* E, unsigned* done, unsigned* yield, unsigned* rn_bank){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	print_histories_kernel <<< blks, NUM_THREADS >>> ( N, isonum, rxn, space, E, done, yield, rn_bank);
	cudaThreadSynchronize();

}

