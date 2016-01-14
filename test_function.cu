#include <cuda.h>
#include <stdio.h>
#include "datadef.h"

__global__ void test_kernel( unsigned N , cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
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
		n_isotopes					= d_xsdata[0].n_isotopes;								
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
	unsigned total_cols = n_isotopes + total_reaction_channels;
	unsigned this_isotope = 2;
	unsigned col_start= n_isotopes + rxn_numbers_total[this_isotope];
	unsigned col_end  = n_isotopes + rxn_numbers_total[this_isotope+1];
	unsigned col = col_start + 10;
	unsigned this_index = row*total_cols+col;

	//
	printf("\n");
	printf("tid %d here isotopes %u this isotope %u\n",tid,n_isotopes,this_isotope);
	printf("energy of grid index %u is %10.8E\n",row,energy_grid[row]);
	printf("column is %u, rxn is %u, total columns %u, index is %u, total xs is %10.8E\n",col,rxn_numbers[col],total_cols,this_index,xs[this_index]);
	printf("Q %6.4E\n",Q[col]);
	printf("awr %6.4E\n",  awr[this_isotope]);
	printf("temp %6.4E\n",temp[this_isotope]);
	printf("scattering dist pointer %p\n",dist_scatter);
	printf("scattering dist pointers, lower %p upper %p\n",dist_scatter[this_index].lower,dist_scatter[this_index].upper);
	if (dist_scatter[this_index].lower != 0x0){
		printf("lower scattering dist, erg %6.8E len %u law %u intt %u\n",dist_scatter[this_index].lower[0].erg,dist_scatter[this_index].lower[0].len,dist_scatter[this_index].lower[0].law,dist_scatter[this_index].lower[0].intt);
		printf("upper scattering dist, erg %6.8E len %u law %u intt %u\n",dist_scatter[this_index].upper[0].erg,dist_scatter[this_index].upper[0].len,dist_scatter[this_index].upper[0].law,dist_scatter[this_index].upper[0].intt);
	}
	else{
		printf("Null dist pointers\n");
	}


}

void test_function( unsigned NUM_THREADS,  unsigned N , cross_section_data* d_xsdata, particle_data* d_particles, tally_data* d_tally){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	test_kernel <<< blks, NUM_THREADS >>> (  N, d_xsdata, d_particles, d_tally );
	cudaThreadSynchronize();

}

