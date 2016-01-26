#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "warp_device.cuh"

__global__ void scatter_multi_kernel(unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

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
		//yield						= d_particles[0].yield;
		weight						= d_particles[0].weight;
		index						= d_particles[0].index;
	}

	// make sure shared loads happen before anything else
	__syncthreads();

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// print and return if wrong
	if ( this_rxn < 11 | this_rxn > 45 ){printf("multiplicity scattering kernel accessing wrong reaction @ dex %u rxn %u\n",tid, this_rxn);return;} 

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in multiplicity scatter!,tid %u rxn %u\n",tid,this_rxn);
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
	float		this_weight		=	weight[  tid];
	//float		this_Q			=	Q[       tid];
	unsigned	rn				=	rn_bank[ tid];
	float		this_awr		=	awr[ this_tope];
	//float		this_temp		=	temp[this_tope];

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
	float  		E_target     		=   0;
	float 		speed_target     	=   sqrtf(2.0*E_target/(this_awr*m_n));
	float  		speed_n          	=   sqrtf(2.0*this_E/m_n);
	float 		E_new				=   0.0;
	float 		sampled_E			=	0.0;
	wfloat3 	v_n_cm, v_t_cm, v_n_lf, v_t_lf, v_cm, hats_new, hats_target, rotation_hat;
	float 		mu, E0, E1, Ek;

	// ensure normalization
	hats_old = hats_old / hats_old.norm2();

	// make speed vectors, assume high enough energy to approximate target as stationary
	v_n_lf = hats_old    * speed_n;
	v_t_lf = hats_target * 0.0;

	// calculate  v_cm
	v_cm = (v_n_lf + (v_t_lf*this_awr))/(1.0+this_awr);

	//transform neutron velocity into CM frame
	v_n_cm = v_n_lf - v_cm;
	v_t_cm = v_t_lf - v_cm;

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

		// sample mu isotropically
		mu  = 2.0*get_rand(&rn)-1.0;

	}
	else if (this_law==44){

		// make sure scatter array is present
		if(dist_scatter == 0x0){
			printf("null pointer, scatter array in continuum !,dex %u rxn %u tope %u E %6.4E \n",this_dex,this_rxn,this_tope,this_E);
			return;
		}

		// compute interpolation factor
		if(f<0){
			printf("DATA NOT WITHIN ENERGY INTERVAL tid %u rxn %u\n",tid,this_rxn);
		}

		// sample tabular on energy, but get index as well as value
		unsigned dist_index = 0;
		E0 = sample_continuous_tablular( 	&dist_index ,
											this_edist.len , 
											this_edist.intt , 
											get_rand(&rn) , 
											this_edist.var , 
											this_edist.pdf, 
											this_edist.cdf );

		//scale it to bins 
		E1 = edist_lower.var[0]                 + f*( edist_upper.var[0]                 - edist_lower.var[0] );
		Ek = edist_lower.var[edist_lower.len-1] + f*( edist_upper.var[edist_upper.len-1] - edist_lower.var[edist_lower.len-1] );
		sampled_E = E1 +(E0-this_edist.var[0])/(this_edist.var[this_edist.len-1]-this_edist.var[0])*(Ek-E1);

		// find correlated mu
		float A 	= this_sdist.var[dist_index];
		float R 	= this_sdist.var[dist_index];
		float rn1 	= get_rand(&rn);
		if( get_rand(&rn)>R ){
			float T = (2.0*rn1-1.0)*sinhf(A);
			mu		= logf(T+sqrtf(T*T+1.0))/A;
		}
		else{
			mu		= logf(rn1*expf(A)+(1.0-rn1)*expf(-A))/A;
		}

	}
	else{

		printf("LAW %u NOT HANDLED IN MULTIPLICITY SCATTER!  rxn %u\n",this_law,this_rxn);

	}

	// rotate direction vector
	hats_old = v_n_cm / v_n_cm.norm2();
	hats_old = hats_old.rotate(mu, get_rand(&rn));

	//  scale to sampled energy
	v_n_cm = hats_old * sqrtf(2.0*sampled_E/m_n);
	
	// transform back to L
	v_n_lf = v_n_cm + v_cm;
	hats_new = v_n_lf / v_n_lf.norm2();
	hats_new = hats_new / hats_new.norm2(); // get higher precision, make SURE vector is length one
	
	// calculate energy in lab frame
	E_new = 0.5 * m_n * v_n_lf.dot(v_n_lf);

	// enforce limits
	if ( E_new <= E_cutoff | E_new > E_max ){
		this_rxn = 998;  // ecutoff code
		printf("c CUTOFF, E = %10.8E\n",E_new);
	}

	// multiply weight by multiplicity
	if(     this_rxn == 22 | this_rxn == 23 | this_rxn == 28 | this_rxn == 29 | this_rxn == 32 | this_rxn == 33 | this_rxn == 34 | this_rxn == 35 | this_rxn == 36 | this_rxn == 44 | this_rxn == 45){
		this_weight = this_weight * 1.0;
	}
	else if(this_rxn == 11 | this_rxn == 16 | this_rxn == 24 | this_rxn == 30 | this_rxn == 41){
		this_weight = this_weight * 2.0;
	}
	else if(this_rxn == 17 | this_rxn == 25 | this_rxn == 42){
		this_weight = this_weight * 3.0;
	}
	else if(this_rxn == 37){
		this_weight = this_weight * 4.0;
	}
	else{
		printf("UNKNOWN MT NUMBER %u IN MULTIPLICITY\n",this_rxn);
	}


	printf("tid %d law %u sampled_E %6.4E mu %6.4E weight %6.4E\n",tid,this_law,sampled_E,mu,this_weight);
	
	// write results
	E[      tid]	=	E_new;
	rn_bank[tid]	=	rn;
	weight[ tid]	=	this_weight;
	space[tid].xhat	=	hats_new.x;
	space[tid].yhat	=	hats_new.y;
	space[tid].zhat	=	hats_new.z;



}

void scatter_multi( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	scatter_multi_kernel <<< blks, NUM_THREADS , 0 , stream >>> ( N, starting_index, d_xsdata, d_particles, d_remap );
	cudaThreadSynchronize();

}
