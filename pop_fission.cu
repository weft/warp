#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void pop_fission_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* scanned){

	// get tid
	int tid = threadIdx.x+blockIdx.x*blockDim.x;

	// declare shared variables					
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	float*				E;		
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			yield;	
	__shared__	unsigned*			index;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){											
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		E							= d_particles[0].E;
		rn_bank						= d_particles[0].rn_bank;
		yield						= d_particles[0].yield;
		index						= d_particles[0].index;
	}

	// load history data
	unsigned		this_dex		=	index[   tid];
	float			this_E			=	E[       tid];
	unsigned		this_yield		=	yield[   tid];
	unsigned		rn				=	rn_bank[ tid];
	unsigned		position		=	scanned[ tid];
	float			this_x			=	space[   tid].x;
	float			this_y			=	space[   tid].y;
	float			this_z			=	space[   tid].z;

	// make sure shared loads happen before anything else (epecially returns)
	__syncthreads();

	// return immediately if out of bounds
	if (tid >= N){return;}

	// check yield
	if (yield[tid]==0){
		return;
	}

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in continuum scatter!,tid %u rxn %u\n",tid);
		return;
	}

	//constants
	const float  	pi			=   3.14159265359;

	// internal kernel variables
	unsigned	data_dex 			=	0;
	unsigned	this_law			=	0;
	float 		sampled_E			=	0.0;
	float 		phi, mu, E0, f;

	// pick upper or lower via stochastic mixing
	dist_data	this_edist, this_sdist;
	dist_data	sdist_lower	=	dist_scatter[this_dex].lower[0];
	dist_data	sdist_upper	=	dist_scatter[this_dex].upper[0];
	dist_data	edist_lower	=	dist_energy[ this_dex].lower[0];
	dist_data	edist_upper	=	dist_energy[ this_dex].upper[0];
	f			=	(this_E - edist_lower.erg) / (edist_upper.erg - edist_lower.erg);

	// write new histories for this yield number
	for(unsigned k=0 ; k < this_yield ; k++ ){

		// do seperate stochastic mixing for each neutron
		if( get_rand(&rn)>f ){
		this_edist	=	edist_lower;
		this_sdist	=	sdist_lower;
		}
		else{
			this_edist = edist_upper;
			this_sdist = sdist_upper;
		}
		this_law	=	this_edist.law;
		
		//get proper data index
		data_dex = position+k;
		
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
			sampled_E = scale_to_bins(	f, E0, 
										 this_edist.var[0],  this_edist.var[ this_edist.len-1], 
										edist_lower.var[0], edist_lower.var[edist_lower.len-1], 
										edist_upper.var[0], edist_upper.var[edist_upper.len-1] );

			// check errors
			if (!isfinite(sampled_E) | sampled_E<=0.0){
				printf("Fission pop mis-sampled tid %i data_dex %u E %6.4E... setting to 2.5\n",tid,data_dex,sampled_E);
				sampled_E = 2.5;
			}
			
			// sample mu/phi isotropically
			mu  = 2.0*get_rand(&rn)-1.0;
			phi = 2.0*pi*get_rand(&rn);

		}
		else{
			printf("LAW %u NOT HANDLED IN FISSION POP!\n",this_law);
		}

		// set data
		E[     data_dex ] 				= sampled_E;
		space[ data_dex ].x				= this_x;
		space[ data_dex ].y				= this_y;
		space[ data_dex ].z				= this_z;
		space[ data_dex ].xhat			= sqrtf(1.0-(mu*mu))*cosf(phi);
		space[ data_dex ].yhat			= sqrtf(1.0-(mu*mu))*sinf(phi); 
		space[ data_dex ].zhat			= mu;
		space[ data_dex ].enforce_BC	= 0;
		space[ data_dex ].surf_dist		= 999999.0;
		
		if(data_dex<=9){printf("array index %u, E = % 6.4E d_fissile_energy[ data_dex ] = % 6.4E\n",data_dex,sampled_E,d_fissile_energy[ data_dex ]);}

	}

	// write current seed out
	rn_bank[tid] = rn;

}

void pop_fission( unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_scanned ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_fission_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_particles, d_scanned);
	check_cuda(cudaThreadSynchronize());

}

