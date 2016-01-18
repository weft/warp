#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void find_E_grid_index_kernel(unsigned N, cross_section_data* d_xsdata , unsigned* remap, float* E , unsigned * index, unsigned* rxn){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	// declare shared variables
	//__shared__ 	unsigned			n_isotopes;				
	__shared__ 	unsigned			energy_grid_len;		
	//__shared__ 	unsigned			total_reaction_channels;
	//__shared__ 	unsigned*			rxn_numbers;			
	//__shared__ 	unsigned*			rxn_numbers_total;		
	__shared__ 	float*				energy_grid;			
	//__shared__ 	float*				Q;						
	//__shared__ 	float*				xs;						
	//__shared__ 	float*				awr;					
	//__shared__ 	float*				temp;					
	//__shared__ 	dist_container*		dist_scatter;			
	//__shared__ 	dist_container*		dist_energy; 

	// have thread 1 copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		//n_isotopes					= d_xsdata[0].n_isotopes;								
		energy_grid_len				= d_xsdata[0].energy_grid_len;				
		//total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		//rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		//rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		energy_grid 				= d_xsdata[0].energy_grid;						
		//Q 							= d_xsdata[0].Q;												
		//xs 							= d_xsdata[0].xs;												
		//awr 						= d_xsdata[0].awr;										
		//temp 						= d_xsdata[0].temp;										
		//dist_scatter 				= d_xsdata[0].dist_scatter;						
		//dist_energy 				= d_xsdata[0].dist_energy;  
	}

	// return if terminated
	unsigned this_rxn=rxn[tid_in];
	if (this_rxn>=900){return;}

	//remap
	int tid=remap[tid_in];

	// load data
	float value = E[tid];
	//unsigned donesearching = 0;
	unsigned cnt  = 0;
	unsigned powtwo = 2;
	int dex  = (energy_grid_len-1) / 2;  //N_energiesgth starts at 1, duh


	for(cnt=0;cnt<=30;cnt++){
		powtwo = powtwo * 2;
		if      ( 	energy_grid[dex]   <= value && 
				energy_grid[dex+1] >  value ) { break; }
		else if ( 	energy_grid[dex]   >  value ) { dex  = dex - ((energy_grid_len / powtwo) + 1) ;}  // +1's are to do a ceiling instead of a floor on integer division
		else if ( 	energy_grid[dex]   <  value ) { dex  = dex + ((energy_grid_len / powtwo) + 1) ;}

		if(cnt==30){
			printf("binary search iteration overflow! %p %d % 10.8f tid=%u\n",energy_grid,energy_grid_len,value,tid);
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


	//write output index
	index[tid]=dex;

}

void find_E_grid_index(unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, unsigned* d_remap, float* d_E , unsigned * d_index , unsigned* d_rxn){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	find_E_grid_index_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_remap,  d_E , d_index , d_rxn);
	cudaThreadSynchronize();

}

