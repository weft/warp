#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void pop_fission_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, spatial_data* d_fissile_points, float* d_fissile_energy, unsigned* scanned){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}

	// declare shared variables					
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			yield;	
	__shared__	unsigned*			index;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){											
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		rn_bank						= d_particles[0].rn_bank;
		yield						= d_particles[0].yield;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	//constants
	const float  	pi			=   3.14159265359;

	// load history data
	spatial_data	this_space		=	space[   tid];
	unsigned		this_rxn		=	rxn[     tid];
	unsigned		this_dex		=	index[   tid];
	float			this_E			=	E[       tid];
	unsigned		this_yield		=	yield[   tid];
	unsigned		rn				=	rn_bank[ tid];
	unsigned		position		=	scanned[ tid];

	// internal kernel variables
	unsigned	data_dex 			=	0;
	unsigned	this_law			=	0;
	float 		sampled_E			=	0.0;
	float 		phi, mu, E0, E1, Ek, f;

	// check yield
	if (yield[tid]==0){
		return;
	}

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in continuum scatter!,tid %u rxn %u\n",tid,this_rxn);
		return;
	}

	// pick upper or lower via stochastic mixing
	dist_data	this_edist, this_sdist;
	dist_data	sdist_lower	=	dist_scatter[this_dex].lower[0];
	dist_data	sdist_upper	=	dist_scatter[this_dex].upper[0];
	dist_data	edist_lower	=	dist_energy[ this_dex].lower[0];
	dist_data	edist_upper	=	dist_energy[ this_dex].upper[0];
	f			=	(this_E - edist_lower.erg) / (edist_upper.erg - edist_lower.erg);
	if( get_rand(&rn)>f ){
		this_edist	=	edist_lower;
		this_sdist	=	sdist_lower;
	}
	else{
		this_edist = edist_upper;
		this_sdist = sdist_upper;
	}
	this_law	=	this_edist.law;

	// write new histories for this yield number
	for(unsigned k=0 ; k < this_yield ; k++ ){
		
		//get proper data index
		data_dex = position+k ;
		
		// sample dist
		if (this_law ==4 ){
	
			// sample continuous tabular
			E0 = sample_continuous_tablular( 	this_edist.len , 
												this_edist.intt , 
												get_rand(&rn) , 
												this_edist.var , 
												this_edist.pdf, 
												this_edist.cdf );
			//scale it to bins 
			E1 = edist_lower.var[0]                 + f*( edist_upper.var[0]                 - edist_lower.var[0] );
			Ek = edist_lower.var[edist_lower.len-1] + f*( edist_upper.var[edist_upper.len-1] - edist_lower.var[edist_lower.len-1] );
			sampled_E = E1 +(E0-this_edist.var[0])/(this_edist.var[this_edist.len-1]-this_edist.var[0])*(Ek-E1);

			// sample mu/phi isotropically
			mu  = 2.0*get_rand(&rn)-1.0;
			phi = 2.0*pi*get_rand(&rn);

		}
		else{
			printf("LAW %u NOT HANDLED IN FISSION POP!  rxn %u\n",this_law,this_rxn);
		}

		// set data
		d_fissile_energy[ data_dex ] 				= sampled_E;
		d_fissile_points[ data_dex ].x 				= this_space.x;
		d_fissile_points[ data_dex ].y 				= this_space.y;
		d_fissile_points[ data_dex ].z 				= this_space.z;
		d_fissile_points[ data_dex ].xhat 			= sqrtf(1.0-(mu*mu))*cosf(phi);
		d_fissile_points[ data_dex ].yhat 			= sqrtf(1.0-(mu*mu))*sinf(phi); 
		d_fissile_points[ data_dex ].zhat 			= mu;
		d_fissile_points[ data_dex ].enforce_BC 	= 0;
		d_fissile_points[ data_dex ].surf_dist 		= 999999.0;

	}

	// write current seed out
	rn_bank[tid] = rn;

}

void pop_fission( unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, spatial_data* d_fissile_points, float* d_fissile_energy, unsigned* scanned){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_fission_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_particles, d_fissile_points, d_fissile_energy, scanned);
	check_cuda(cudaThreadSynchronize());

}

