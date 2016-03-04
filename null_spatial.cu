#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"
#include "wfloat3.h"

__global__ void null_spatial_kernel(unsigned N, spatial_data* space){   

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;} 

	space[k].x				= 0.0;
	space[k].y				= 0.0;
	space[k].z				= 0.0;
	space[k].xhat			= 0.0;
	space[k].yhat			= 0.0;
	space[k].zhat			= 0.0;
	space[k].surf_dist		= 10000.0;
	space[k].enforce_BC		= 0.0;
	space[k].norm[0]		= 1.0;
	space[k].norm[1]		= 0.0;
	space[k].norm[2]		= 0.0;

}

void null_spatial(unsigned NUM_THREADS, unsigned N, spatial_data* space){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	null_spatial_kernel <<< blks, NUM_THREADS >>> (  N, space );
	check_cuda(cudaThreadSynchronize());

}

