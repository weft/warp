#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void find_E_grid_index_kernel(unsigned N, unsigned N_energies , unsigned* remap, float * main_E_grid, float* E , unsigned * index, unsigned* rxn){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>801){return;}

	//remap
	int tid=remap[tid_in];

	// load data
	float value = E[tid];
	//unsigned donesearching = 0;
	unsigned cnt  = 0;
	unsigned powtwo = 2;
	int dex  = (N_energies-1) / 2;  //N_energiesgth starts at 1, duh


	for(cnt=0;cnt<=30;cnt++){
		powtwo = powtwo * 2;
		if      ( 	main_E_grid[dex]   <= value && 
				main_E_grid[dex+1] >  value ) { break; }
		else if ( 	main_E_grid[dex]   >  value ) { dex  = dex - ((N_energies / powtwo) + 1) ;}  // +1's are to do a ceiling instead of a floor on integer division
		else if ( 	main_E_grid[dex]   <  value ) { dex  = dex + ((N_energies / powtwo) + 1) ;}

		if(cnt==30){
			printf("binary search iteration overflow! %p %d % 10.8f tid=%u\n",main_E_grid,N_energies,value,tid);
			dex=0;
		}

		// edge checks... fix later???
		if(dex<0){
			dex=0;
		}
		if(dex>=N_energies){
			dex=N_energies-1;
		}
	}


	//write output index
	index[tid]=dex;

}

void find_E_grid_index(unsigned NUM_THREADS, unsigned N, unsigned N_energies, unsigned* active, float * main_E_grid, float* E , unsigned * index , unsigned* rxn){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	find_E_grid_index_kernel <<< blks, NUM_THREADS >>> ( N, N_energies, active, main_E_grid,  E , index , rxn);
	cudaThreadSynchronize();

}

