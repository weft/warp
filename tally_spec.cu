#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void tally_spec_kernel(unsigned N, unsigned Ntally, unsigned tally_cell,  unsigned* remap, source_point* space, float* E, float * tally_score, float * tally_square, unsigned * tally_count, unsigned* done, unsigned* cellnum, unsigned* rxn){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// remap to active
	int tid = remap[tid_in];
	unsigned this_rxn = rxn[tid_in];   // rxn is sorted along with remap vector, must use tid_in.

	// return if not cell or for some reason marked as done
	if (cellnum[tid]!=tally_cell){return;}
	if (this_rxn>=900 | this_rxn==800){return;}

	//int k;
	float 		my_E   			= E[tid];
	float 		macro_t 		= space[tid].macro_t;
	unsigned 	my_bin_index 	= 0;

	const float Emax 	= 20.00000;
	const float Emin 	=  1.0e-11;

	// determine bin number
	my_bin_index = logf(my_E/Emin)/logf(Emax/Emin)*(Ntally);

	//score the bins atomicly, could be bad if many neutrons are in a single bin since this will serialize their operations
	atomicAdd(&tally_score [my_bin_index], 1.0/macro_t );
	atomicAdd(&tally_square[my_bin_index], 1.0/(macro_t * macro_t));
	atomicInc(&tally_count [my_bin_index], 4294967295);

	//printf("%6.4E\n",macro_t);


}

void tally_spec(unsigned NUM_THREADS,  unsigned N, unsigned Ntally, unsigned tally_cell, unsigned* remap, source_point * space, float* E, float * tally_score, float * tally_square, unsigned * tally_count, unsigned* done, unsigned* cellnum, unsigned* rxn){
	
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	tally_spec_kernel <<< blks, NUM_THREADS >>> ( N, Ntally, tally_cell, remap, space, E, tally_score, tally_square, tally_count, done, cellnum, rxn);
	cudaThreadSynchronize();

}

