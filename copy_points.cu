#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void copy_points_kernel( unsigned Nout, unsigned * Nvalid , unsigned current_index , unsigned * to_valid, source_point * positions_out , source_point * positions_in, float*E_out, float*E_in  ){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= Nvalid[0]){return;}

	unsigned index_in  = to_valid[tid];
	unsigned index_out = current_index + tid;
	if (index_out>=Nout){index_out=index_out-Nout;} //wrap to start
	//printf("index out = %u\n",index_out);

	// copy points
	positions_out[index_out].x    	= positions_in[index_in].x; 
	positions_out[index_out].y    	= positions_in[index_in].y; 
	positions_out[index_out].z    	= positions_in[index_in].z; 
	positions_out[index_out].xhat 	= positions_in[index_in].xhat; 
	positions_out[index_out].yhat 	= positions_in[index_in].yhat; 
	positions_out[index_out].zhat 	= positions_in[index_in].zhat;
	E_out[index_out] 				= E_in[index_in];

	//printf("good point %6.4E %6.4E %6.4E\n",positions_out[index_out].x,positions_out[index_out].y,positions_out[index_out].z);


}

void copy_points( unsigned NUM_THREADS,  unsigned Nout , unsigned * Nvalid,  unsigned current_index , unsigned * to_valid , source_point * positions_out , source_point * positions_in, float*E_out, float*E_in){

	unsigned blks = ( Nout + NUM_THREADS - 1 ) / NUM_THREADS;

	copy_points_kernel <<< blks, NUM_THREADS >>> (  Nout , Nvalid,  current_index , to_valid , positions_out , positions_in , E_out, E_in);
	cudaThreadSynchronize();

}

