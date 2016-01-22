#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "warp_device.cuh"

__global__ void fission_kernel(unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}       

	// declare shared variables
	//__shared__ 	unsigned			n_isotopes;				
	//__shared__ 	unsigned			energy_grid_len;		
	//__shared__ 	unsigned			total_reaction_channels;
	//__shared__ 	unsigned*			rxn_numbers;			
	//__shared__ 	unsigned*			rxn_numbers_total;		
	//__shared__ 	float*				energy_grid;			
	//__shared__ 	float*				rxn_Q;						
	//__shared__ 	float*				xs;						
	__shared__ 	float*				awr;					
	//__shared__ 	float*				temp;					
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	//__shared__	float*				Q;		
	__shared__	unsigned*			rn_bank;
	//__shared__	unsigned*			cellnum;
	//__shared__	unsigned*			matnum;	
	__shared__	unsigned*			isonum;	
	//__shared__	unsigned*			yield;	
	//__shared__	float*				weight;	
	__shared__	unsigned*			index;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		//n_isotopes					= d_xsdata[0].n_isotopes;								
		//energy_grid_len				= d_xsdata[0].energy_grid_len;				
		//total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		//rxn_numbers 				= d_xsdata[0].rxn_numbers;						
		//rxn_numbers_total			= d_xsdata[0].rxn_numbers_total;					
		//energy_grid 				= d_xsdata[0].energy_grid;						
		//rxn_Q 						= d_xsdata[0].Q;												
		//xs 							= d_xsdata[0].xs;												
		awr 						= d_xsdata[0].awr;										
		//temp 						= d_xsdata[0].temp;										
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		//Q							= d_particles[0].Q;	
		rn_bank						= d_particles[0].rn_bank;
		//cellnum						= d_particles[0].cellnum;
		//matnum						= d_particles[0].matnum;
		isonum						= d_particles[0].isonum;
		//yield						= d_particles[0].yield;
		//weight						= d_particles[0].weight;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// print and return if wrong
	if ( this_rxn!=91 ){printf("level scattering kernel accessing wrong reaction @ dex %u rxn %u\n",tid, this_rxn);return;} 

	//constants
	//const float  	pi			=   3.14159265359;
	const float  	m_n			=   1.00866491600;		// u
	const float  	E_cutoff	=   1e-11;				// MeV
	const float  	E_max		=   20.0;				// MeV
	//const float		kb			=	8.617332478e-11;	// MeV/k

	// load history data
	wfloat3		hats_old(space[tid].xhat,space[tid].yhat,space[tid].zhat);
	unsigned	this_tope		=	isonum[  tid];
	unsigned	this_dex		=	index[   tid];
	float		this_E			=	E[       tid];
	//float		this_Q			=	Q[       tid];
	unsigned	rn				=	rn_bank[ tid];
	float		this_awr		=	awr[ this_tope];
	//float		this_temp		=	temp[this_tope];

	

}

void fission( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	fission_kernel <<< blks, NUM_THREADS , 0 , stream >>> ( N, starting_index, d_xsdata, d_particles, d_remap );
	cudaThreadSynchronize();

}

