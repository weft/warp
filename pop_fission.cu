#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void pop_fission_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_scanned){

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

	// make sure shared loads happen before anything else (epecially returns)
	__syncthreads();

	// load history data
	unsigned		this_dex		=	index[    tid];
	float			this_E			=	E[        tid];
	unsigned		this_yield		=	yield[    tid];
	unsigned		rn				=	rn_bank[  tid];
	float			this_x			=	space[    tid].x;
	float			this_y			=	space[    tid].y;
	float			this_z			=	space[    tid].z;

	// get array position from prefix scan
	unsigned	position	=	d_scanned[tid];

	// make sure individual loads happen before anything else?
	__syncthreads();

	// return immediately if out of bounds
	if (tid >= N){return;}

	// check yield
	if (this_yield==0){
		return;
	}

	// another yield check
	if((d_scanned[tid+1]-d_scanned[tid]) == 0){
		printf("NOT RIGHT! \n");
		return;
	}

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in continuum scatter!,tid %u\n",tid);
		return;
	}

	//constants
	const float  	pi			=   3.14159265359;

	// internal kernel variables
	float		nu_t0				=	0.0;
	float		nu_t1				=	0.0;
	float		nu_d0				=	0.0;
	float		nu_d1				=	0.0;
	float		beta				=	0.0;
	float		e0					=	0.0;
	float		e1					=	0.0;
	unsigned	data_dex 			=	0;
	float 		sampled_E			=	0.0;
	float 		phi, mu, E0, f, rn1;
	unsigned	this_law, this_len, this_intt, upper_len, lower_len, pre_index, pre_position;
	float		*this_var, *this_cdf, *this_pdf, *upper_var, *lower_var;

	// load dist info
	dist_data	this_edist, this_sdist;
	dist_data	sdist_lower	=	dist_scatter[this_dex].lower[0];
	dist_data	sdist_upper	=	dist_scatter[this_dex].upper[0];
	dist_data	edist_lower	=	dist_energy[ this_dex].lower[0];
	dist_data	edist_upper	=	dist_energy[ this_dex].upper[0];

	// copy nu values, energy points from dist, t is len, d is law
	memcpy(&nu_t0	, &sdist_lower.len, 1*sizeof(float));
	memcpy(&nu_t1	, &sdist_upper.len, 1*sizeof(float));
	memcpy(&nu_d0	, &sdist_lower.law, 1*sizeof(float));
	memcpy(&nu_d1	, &sdist_upper.law, 1*sizeof(float));
	memcpy(&e0		, &sdist_lower.erg, 1*sizeof(float));
	memcpy(&e1		, &sdist_upper.erg, 1*sizeof(float));

	// get interpolated beta value, beta = nu_d / nu_t
	beta	=	interpolate_linear_energy( this_E, e0, e1, nu_d0, nu_d1 ) / 
				interpolate_linear_energy( this_E, e0, e1, nu_t0, nu_t1 );

	// write new histories for this yield number
	for(unsigned k=0 ; k < this_yield ; k++ ){

		//get proper data index
		data_dex = position+k;

		// check if this neutron is delayed or prompt
		if ( get_rand(&rn) > beta ){

			// do individual stochastic mixing for this prompt neutron
			// pick upper or lower edist via stochastic mixing
			f	=	(this_E - edist_lower.erg) / (edist_upper.erg - edist_lower.erg);
			if( get_rand(&rn) > f ){
				this_edist	=	edist_lower;
			}
			else{
				this_edist	=	edist_upper;
			}

			// set pointers and parameters
			this_law	=	this_edist.law;
			this_len 	=	this_edist.len;
			this_intt	=	this_edist.intt;
			this_var	=	this_edist.var;
			this_cdf	=	this_edist.cdf;
			this_pdf	=	this_edist.pdf;
			upper_var	=	edist_upper.var;
			lower_var	=	edist_lower.var;
			upper_len	=	edist_upper.len;
			lower_len	=	edist_lower.len;

		}
		else{

			// pick upper or lower sdist (contains the delayed data) via stochastic mixing
			f	=	0.0;//(this_E - sdist_lower.var[0]) / (sdist_upper.erg - sdist_lower.erg);
			if( get_rand(&rn) > f ){
				this_sdist	=	sdist_lower;
			}
			else{
				this_sdist	=	sdist_upper;
			}

			// decode precursor intt, 100 place
			this_intt	=	(this_sdist.intt%1000-this_sdist.intt%100)/100;

			// decode precursor law, 1000 place
			this_law	=	 (this_sdist.intt%10000-this_sdist.intt%1000)/1000;

			// sample which precursor neutron is from
			rn1 = get_rand(&rn);
			for( pre_index=0; pre_index<6; pre_index++ ){ 
				if ( rn1 <= this_sdist.var[pre_index+1] ){
					break;
				}
			}

			// get position of data in vector and vector length
			pre_position	=	(unsigned) this_sdist.pdf[pre_index];  									// haha preposition...
			this_len		=	(unsigned) this_sdist.pdf[pre_index+1] - (unsigned) this_sdist.pdf[pre_index];

			// get pointers to sampled data
			this_var	=	&this_sdist.cdf[pre_position];
			this_cdf	=	&this_sdist.cdf[pre_position + ((unsigned)this_sdist.pdf[6])   ];   // last value in cdf if the total length of the combined 6-vectors
			this_pdf	=	&this_sdist.cdf[pre_position + ((unsigned)this_sdist.pdf[6])*2 ];
			upper_var	=	&this_sdist.cdf[pre_position];
			lower_var	=	&this_sdist.cdf[pre_position];
			upper_len	=	this_len;
			lower_len	=	this_len;

			printf("DELAYED this_E %6.4E f %6.4E  pre_index %u pre_position %u this_len %u this_var[0] %6.4E this_cdf[0] %6.4E this_pdf[0] %6.4E\n",this_E,f,pre_index,pre_position,this_len,this_var[0],this_cdf[0],this_pdf[0]);

		}
			
		// sample dist, passing the parameters/pointers to the sampled delayed/prompt emission data
		if (this_law ==4 ){
	
			// sample continuous tabular
			E0 = sample_continuous_tablular( 	this_len , 
												this_intt , 
												get_rand(&rn) , 
												this_var , 
												this_cdf, 
												this_pdf );
			//scale it to bins 
			sampled_E = scale_to_bins(	f, E0, 
										 this_var[0],  this_var[ this_len-1], 
										lower_var[0], lower_var[lower_len-1], 
										upper_var[0], upper_var[upper_len-1] );
		
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
		
		//if(data_dex<=9){printf("array index %u, E = % 6.4E d_fissile_energy[ data_dex ] = % 6.4E\n",data_dex,sampled_E,E[ data_dex ]);}

	}

	// write current seed out
	rn_bank[tid] = rn;

}

void pop_fission( unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_scanned ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_fission_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_particles, d_scanned);
	check_cuda(cudaThreadSynchronize());

}

