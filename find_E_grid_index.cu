#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void find_E_grid_index_kernel(unsigned N, cross_section_data* d_xsdata , unsigned* remap, float* E , unsigned * index, unsigned* rxn){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// declare shared variables			
	__shared__ 	unsigned			energy_grid_len;				
	__shared__ 	float*				energy_grid;			

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		energy_grid_len				= d_xsdata[0].energy_grid_len;								
		energy_grid 				= d_xsdata[0].energy_grid;						
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>=900){return;}

	// remap
	int tid=remap[tid_in];

	// load data
	float value = E[tid];

	// init local
	unsigned cnt  = 0;
	unsigned powtwo = 2;
	int dex  = (energy_grid_len-1) / 2;  // N_energies starts at 1, duh

	// check edges
	float	grid_max = energy_grid[energy_grid_len-1];
	float	grid_min = energy_grid[0];
	if (value < grid_min){
		dex=-2;   // encoding for under min interpolation (2^32-2, or -2)
	}
	else if (value > grid_max){
		dex=-1;   // encoding for over max interpolation (2^32-1, or -1)
	}
	else{
		// do the search
		for(cnt=0;cnt<=30;cnt++){
			powtwo = powtwo * 2;
			if      ( 	energy_grid[dex]   <= value && 
					energy_grid[dex+1] >  value ) { break; }
			else if ( 	energy_grid[dex]   >  value ) { dex  = dex - ((energy_grid_len / powtwo) + 1) ;}  // +1's are to do a ceiling instead of a floor on integer division
			else if ( 	energy_grid[dex]   <  value ) { dex  = dex + ((energy_grid_len / powtwo) + 1) ;}
	
			if(cnt==30){
				printf("binary search iteration overflow! %p len %d val % 6.4E tid=%u rxn=%u\n",energy_grid,energy_grid_len,value,tid,this_rxn);
				dex=0;
			}
	
			// edge checks... fix later???
			if(dex<0){
				dex=0;
			}
			if(dex>=energy_grid_len){
				dex=energy_grid_len-1;
			}
		}
	}

	//write output index
	index[tid]=dex;

	//printf("rxn %u, remap[%i]=%u E %6.4E dex %u\n",this_rxn,tid_in,tid, value, dex);

}


/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS    - the number of threads to run per thread block
 * @param[in]    N              - the total number of threads to launch on the grid
 * @param[in]    d_xsdata       - device pointer to cross section data pointer array 
 * @param[in]    d_remap        - device pointer to data remapping vector
 * @param[in]    d_E            - device pointer to energy data array
 * @param[in]    d_index        - device pointer to index array (stores the unionized grid index of the current energy)
 * @param[in]    d_rxn          - device pointer of the reaction number array
 */ 
void find_E_grid_index(unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, unsigned* d_remap, float* d_E , unsigned * d_index , unsigned* d_rxn){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	find_E_grid_index_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_remap,  d_E , d_index , d_rxn);
	cudaThreadSynchronize();

}

