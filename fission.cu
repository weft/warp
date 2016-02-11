#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void fission_kernel(unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
	__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	float*				energy_grid;
	__shared__ 	dist_container*		dist_scatter;			
	__shared__	unsigned*			rxn;	
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			yield;	
	__shared__	float*				weight;
	__shared__	float*				E;	
	__shared__	unsigned*			index;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;								
		energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		energy_grid 				= d_xsdata[0].energy_grid;
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		rxn							= d_particles[0].rxn;
		rn_bank						= d_particles[0].rn_bank;
		yield						= d_particles[0].yield;
		weight						= d_particles[0].weight;
		index						= d_particles[0].index;
		E							= d_particles[0].E;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// print and return if wrong
	if ( this_rxn != 818 & this_rxn != 819 & this_rxn != 820 & this_rxn != 821 & this_rxn != 838 ){printf("fission kernel accessing wrong reaction @ dex %u dex_in %u rxn %u\n",tid, tid_in,this_rxn);return;} 

	// load history data
	unsigned	this_dex		=	index[  tid];
	unsigned	rn				=	rn_bank[tid];
	float		this_weight		=	weight[ tid];
	float		this_E			=	E[      tid];

	// local variables, load nu from scattering dist variables
	if (dist_scatter[this_dex].lower==0x0){
		printf("scatter pointer for rxn %d is null!\n",this_rxn);
	}
	float		nu0			=	0.0;
	float		nu1			=	0.0;
	float		nu			=	0.0;
	unsigned	inu			=	0;
	unsigned	this_yield	=	0;
	unsigned 	n_columns	= n_isotopes + total_reaction_channels;

	// copy nu from array, lower erg is nu_t, upper erg is nu_p
	memcpy(&nu0, &dist_scatter[this_dex          ].lower, 1*sizeof(float));
	memcpy(&nu1, &dist_scatter[this_dex+n_columns].lower, 1*sizeof(float));

	// interpolate nu
	float	f	=	(this_E - energy_grid[this_dex]) / (energy_grid[this_dex+1] - energy_grid[this_dex]);
	nu	=	f*(nu1 - nu0) + nu0;

	// check nu
	if (nu==0.0){
		nu=2.8;
		printf("something is wrong with fission yields, nu = %6.4E, guessing %4.2f, rxn %u\n",0.0,nu,this_rxn); 
	}

	//  multiply nu by weight
	nu = this_weight * nu;
	if (this_weight<1.0){
		printf("weight = 0!\n");
	}
	//printf("nu %6.4E\n",nu);

	// get integer part
	inu = (unsigned) nu;
	
	// sample floor or ceil based on fractional part
	if((float)inu+get_rand(&rn) <= nu){
		this_yield = inu+1;
	}
	else{
		this_yield = inu;
	}

	// put in 900 block to terminate on next sort
 	this_rxn += 100;

	//printf("tid %d rxn %u wgt %6.4E yield %u\n", tid, this_rxn, this_weight, this_yield);

	// write 
	yield[  tid]					=	this_yield;
	rn_bank[tid]					=	rn;  
	rxn[starting_index + tid_in]	=	this_rxn;

}

void fission( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	fission_kernel <<< blks, NUM_THREADS , 0 , stream >>> ( N, starting_index, d_xsdata, d_particles, d_remap );
	check_cuda(cudaThreadSynchronize());

}

