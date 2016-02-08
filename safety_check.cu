#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"
#include "wfloat3.h"

__global__ void safety_check_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){   

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
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	__shared__	float*				Q;		
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			cellnum;
	__shared__	unsigned*			matnum;	
	__shared__	unsigned*			isonum;	
	__shared__	unsigned*			yield;	
	__shared__	float*				weight;	
	__shared__	unsigned*			index;	

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
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		Q							= d_particles[0].Q;	
		rn_bank						= d_particles[0].rn_bank;
		cellnum						= d_particles[0].cellnum;
		matnum						= d_particles[0].matnum;
		isonum						= d_particles[0].isonum;
		yield						= d_particles[0].yield;
		weight						= d_particles[0].weight;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;} 

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// check energy
	float 	this_E = E[tid];
	if (!isfinite(this_E) | this_E < 0.0){
		printf("INVALID ENERGY, tid %u tid_in %u rxn %u, E % 6.4E\n",tid,tid_in,this_rxn,this_E);
	}

	// check directions
	wfloat3		hats(space[tid].xhat,space[tid].yhat,space[tid].zhat);
	if (!isfinite(hats.x+hats.y+hats.z)){
		printf("INVALID DIRECTIONS, tid %u tid_in %u rxn %u, xhat % 6.4E yhat % 6.4E zhat % 6.4E\n",tid,tid_in,this_rxn,hats.x,hats.y,hats.z);
	}

	// check position
	wfloat3		pos(space[tid].x,space[tid].y,space[tid].z);
	if (!isfinite(pos.x+pos.y+pos.z)){
		printf("INVALID POSITIONS, tid %u tid_in %u rxn %u, x % 6.4E y % 6.4E z % 6.4E\n",tid,tid_in,this_rxn,pos.x,pos.y,pos.z);
	}

}

void safety_check(unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	safety_check_kernel <<< blks, NUM_THREADS >>> (  N, d_xsdata, d_particles, d_remap);
	check_cuda(cudaThreadSynchronize());

}

