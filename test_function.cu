#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void test_kernel( unsigned N , cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// declare shared variables
	__shared__ 	unsigned			isotopes;				
	__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	unsigned*			rxn_numbers;			
	__shared__ 	unsigned*			rxn_numbers_total;		
	__shared__ 	float*				energy_grid;			
	__shared__ 	float*				Q;						
	__shared__ 	float*				xs;						
	__shared__ 	float*				awr;					
	__shared__ 	float*				temp;					
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 

	// have thread 1 copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		isotopes					= d_xsdata[0].isotopes;								
		energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		energy_grid 				= d_xsdata[0].energy_grid;						
		Q 							= d_xsdata[0].Q;												
		xs 							= d_xsdata[0].xs;												
		awr 						= d_xsdata[0].awr;										
		temp 						= d_xsdata[0].temp;										
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		dist_energy 				= d_xsdata[0].dist_energy;  
	}

	// go about your thready business
	unsigned row = energy_grid_len*0.99;
	unsigned total_cols = isotopes + total_reaction_channels;
	unsigned this_isotope = 0;
	unsigned col_start= isotopes + rxn_numbers_total[this_isotope];
	unsigned col_end  = isotopes + rxn_numbers_total[this_isotope+1];

	unsigned col = col_start + 0;

	//
	printf("tid %d here isotopes %u this isotope %u\n",tid,isotopes,this_isotope);
	printf("energy of grid index %u is %10.8E\n",row,energy_grid[row]);
	printf("column is %u, rxn is %u, total columns %u, index is %u, total xs is %10.8E\n",col,rxn_numbers[col],total_cols,row*total_cols+col,xs[row*total_cols+col]);

	printf("scattering dist pointer %p\n",dist_scatter);


}

void test_function( unsigned NUM_THREADS,  unsigned N , cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	test_kernel <<< blks, NUM_THREADS >>> (  N, d_xsdata, d_particles, d_tally );
	cudaThreadSynchronize();

}

