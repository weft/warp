#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"

__global__ void pop_fission_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* scanned){

	int tid = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid >= N){return;}
	if (yield[tid]==0){return;}

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
	__shared__	unsigned*			yield;	
	__shared__	float*				weight;	
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
		yield						= d_particles[0].yield;
		weight						= d_particles[0].weight;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in continuum scatter!,tid %u rxn %u\n",tid,this_rxn);
		return;
	}

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
	unsigned	this_yield		=	yield[   tid];
	float		this_weight		=	weight[  tid];
	unsigned	rn				=	rn_bank[ tid];
	unsigned	position		=	scanned[ tid];
	float		this_awr		=	awr[ this_tope];

	// pick upper or lower via stochastic mixing
	dist_data	this_edist, this_sdist;
	dist_data	sdist_lower	=	dist_scatter[this_dex].lower[0];
	dist_data	sdist_upper	=	dist_scatter[this_dex].upper[0];
	dist_data	edist_lower	=	dist_energy[ this_dex].lower[0];
	dist_data	edist_upper	=	dist_energy[ this_dex].upper[0];
	unsigned	this_law;
	float		f			=	(this_E - edist_lower.erg) / (edist_upper.erg - edist_lower.erg);
	if( get_rand(&rn)>f ){
		this_edist	=	edist_lower;
		this_sdist	=	sdist_lower;
	}
	else{
		this_edist = edist_upper;
		this_sdist = sdist_upper;
	}
	this_law	=	this_edist.law;

	// internal kernel variables
	float 		E_new				=   0.0;
	float 		sampled_E			=	0.0;
	float 		mu, E0, E1, Ek;

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

	//check limits
	if (sampled_E >= Emax){sampled_E = Emax * 0.9;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
	if (sampled_E <= Emin){sampled_E = Emin * 1.1;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}

	// set self data
	E    [ tid ] 		= sampled_E;
	space[ tid ].xhat 	= sqrtf(1.0-(mu*mu))*cosf(phi);
	space[ tid ].yhat 	= sqrtf(1.0-(mu*mu))*sinf(phi); 
	space[ tid ].zhat 	= mu;
	done [ tid ] 		= 0;
	yield[ tid ] 		= 0;
	rxn  [ tid ] 		= 0;//this_rxn;
	

	for(k=0 ; k < this_yield-1 ; k++ ){
		//get proper data index
		data_dex=completed[position+k];

		//make sure data is done
		if(!done[data_dex]){printf("overwriting into active data!\n");}
		
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
	
		//check limits
		if (sampled_E >= Emax){sampled_E = Emax * 0.9;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
		if (sampled_E <= Emin){sampled_E = Emin * 1.1;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}

		// set data
		E    [ data_dex ] 				= sampled_E;
		space[ data_dex ].x 			= this_space.x;
		space[ data_dex ].y 			= this_space.y;
		space[ data_dex ].z 			= this_space.z;
		space[ data_dex ].xhat 			= sqrtf(1.0-(mu*mu))*cosf(phi);
		space[ data_dex ].yhat 			= sqrtf(1.0-(mu*mu))*sinf(phi); 
		space[ data_dex ].zhat 			= mu;
		space[ data_dex ].enforce_BC 	= 0;
		space[ data_dex ].surf_dist 	= 999999.0;
		space[ data_dex ].macro_t 		= 9999999.0;
		done [ data_dex ] 				= 0;
		yield[ data_dex ] 				= 0;
		rxn  [ data_dex ]				= 0;

	}

	// write current seed out
	rn_bank[tid] = rn;

}

void pop_fission( unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* scanned){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_fission_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_particles, scanned);
	cudaThreadSynchronize();

}

