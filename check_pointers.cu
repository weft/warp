#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"
#include "wfloat3.h"

__global__ void check_pointers_kernel(unsigned N, unsigned dex0, unsigned dex1, cross_section_data* d_xsdata ){   

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
	__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	unsigned*			rxn_numbers;			
	__shared__ 	unsigned*			rxn_numbers_total;		
	__shared__ 	float*				energy_grid;			
	__shared__ 	float*				rxn_Q;						
	__shared__ 	float*				xs;						
	__shared__ 	float*				awr;					
	__shared__ 	float*				temp;					
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;								
		energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		energy_grid 				= d_xsdata[0].energy_grid;						
		rxn_Q 						= d_xsdata[0].Q;												
		xs 							= d_xsdata[0].xs;												
		awr 						= d_xsdata[0].awr;										
		temp 						= d_xsdata[0].temp;										
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		dist_energy 				= d_xsdata[0].dist_energy; 
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return immediately if out of bounds
	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;} 

	printf("INDEX %u -> energy  pointer = %p -> lower %p upper %p\n",dex0+tid,dist_energy, dist_energy[ dex0+tid].upper,dist_energy[ dex0+tid].lower);
	printf("INDEX %u -> scatter pointer = %p -> lower %p upper %p\n",dex0+tid,dist_scatter,dist_scatter[dex0+tid].upper,dist_scatter[dex0+tid].lower);

}

void check_pointers(unsigned NUM_THREADS, unsigned dex0, unsigned dex1, cross_section_data* d_xsdata){

	int N = dex1-dex0+1;
	if (N<1){printf("Negative range in check_pointers! dex0 %u dex1 %u -> N = %d\n",dex0,dex1,N);return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	check_pointers_kernel <<< blks, NUM_THREADS >>> ( N, dex0, dex1, d_xsdata );
	check_cuda(cudaThreadSynchronize());

}

