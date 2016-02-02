#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void scatter_level_kernel(unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

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
	__shared__ 	float*				temp;					
	__shared__ 	dist_container*		dist_scatter;			
	//__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	unsigned*			rxn;	
	__shared__	float*				E;		
	__shared__	float*				Q;		
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
		temp 						= d_xsdata[0].temp;										
		dist_scatter 				= d_xsdata[0].dist_scatter;						
		//dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		rxn							= d_particles[0].rxn;
		E							= d_particles[0].E;
		Q							= d_particles[0].Q;	
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
	if ( this_rxn < 50 | this_rxn > 90 ){printf("level scattering kernel accessing wrong reaction @ dex %u dex_in %u rxn %u\n",tid, tid_in, this_rxn);return;} 

	//constants
	const float  	pi			=   3.14159265359;
	const float  	m_n			=   1.00866491600;		// u
	const float		kb			=	8.617332478e-11;	// MeV/k

	// load history data
	wfloat3		hats_old(space[tid].xhat,space[tid].yhat,space[tid].zhat);
	unsigned	this_tope		=	isonum[  tid];
	unsigned	this_dex		=	index[   tid];
	float		this_E			=	E[       tid];
	float		this_Q			=	Q[       tid];
	unsigned	rn				=	rn_bank[ tid];
	float		this_awr		=	awr[ this_tope];
	float		this_temp		=	temp[this_tope];

	// internal kernel variables
	float 		speed_target;     	
	float  		speed_n          	=   sqrtf(2.0*this_E/m_n);
	float 		E_new				=   0.0;
	wfloat3 	v_n_cm, v_t_cm, v_n_lf, v_t_lf, v_cm, hats_new, hats_target, rotation_hat;
	float 		mu, phi, arg;

	// ensure normalization
	hats_old = hats_old / hats_old.norm2();

	// make target isotropic
	mu				=	(2.0*   get_rand(&rn)) - 1.0;
	phi				=	 2.0*pi*get_rand(&rn);
	hats_target.x	=	sqrtf(1.0-(mu*mu))*cosf(phi);
	hats_target.y	=	sqrtf(1.0-(mu*mu))*sinf(phi); 
	hats_target.z	=	mu;

	//sample therm dist if low E
	if(this_E <= 600*kb*this_temp ){
		sample_therm(&rn,&mu,&speed_target,this_temp,this_E,this_awr);
		hats_target = hats_old.rotate(mu, get_rand(&rn));
	}
	else{
		speed_target = 0.0;
	}
	__syncthreads();
	
	// make speed vectors
	v_n_lf = hats_old    * speed_n;
	v_t_lf = hats_target * speed_target;

	// calculate  v_cm
	v_cm = (v_n_lf + (v_t_lf*this_awr))/(1.0+this_awr);

	//transform neutron velocity into CM frame
	v_n_cm = v_n_lf - v_cm;
	v_t_cm = v_t_lf - v_cm;
	
	// sample new azimuthal phi, always isotropic
	phi = 2.0*pi*get_rand(&rn);

	// sample polar cosine mu
	if(dist_scatter == 0x0){
		// isotropic scatter if null
		mu= 2.0*get_rand(&rn)-1.0;
	}
	else{  
		// pick upper or lower via stochastic mixing
		dist_data	this_dist;
		dist_data	dist_lower	=	dist_scatter[this_dex].lower[0];
		dist_data	dist_upper	=	dist_scatter[this_dex].upper[0];
		unsigned	this_law	=	0;
		float		f			=	(this_E - dist_lower.erg) / (dist_upper.erg - dist_lower.erg);
		if( get_rand(&rn)>f ){
			this_dist = dist_lower;
		}
		else{
			this_dist = dist_upper;
		}

		// sample the distribution
		this_law = this_dist.law;
		if (this_law == 3 ){
			mu = sample_continuous_tablular( 	this_dist.len , 
												this_dist.intt , 
												get_rand(&rn) , 
												this_dist.var , 
												this_dist.pdf, 
												this_dist.cdf );
		}
		else{
			printf("law %u not yet implemented in level scttering!\n",this_law);
		}
	}

	// pre rotation directions
	hats_old = v_n_cm / v_n_cm.norm2();
	hats_old = hats_old.rotate(mu, get_rand(&rn));

	// check arg to make sure not negative
	arg = v_n_cm.dot(v_n_cm) + 2.0*this_awr*this_Q/((this_awr+1.0)*m_n);
	if(arg < 0.0) { 
		arg=0.0;
	}
	v_n_cm = hats_old * sqrtf( arg );

	// transform back to L frame
	v_n_lf = v_n_cm + v_cm;
	hats_new = v_n_lf / v_n_lf.norm2();
	hats_new = hats_new / hats_new.norm2();  // get higher precision, make SURE vector is length one
	
	// calculate energy in L frame
	E_new = 0.5 * m_n * v_n_lf.dot(v_n_lf);

	//printf("tid %d E_new %6.4E xhat %6.4E yhat %6.4E zhat %6.4E\n",tid,this_E,hats_new.x,hats_new.y,hats_new.z);

	// write results
	E[      tid]		= E_new;
	space[  tid].xhat	= hats_new.x;
	space[  tid].yhat	= hats_new.y;
	space[  tid].zhat	= hats_new.z;
	rn_bank[tid]		= rn;

}

void scatter_level( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	scatter_level_kernel <<< blks, NUM_THREADS , 0 , stream >>> (  N, starting_index, d_xsdata, d_particles, d_remap);
	check_cuda(cudaThreadSynchronize());

}

