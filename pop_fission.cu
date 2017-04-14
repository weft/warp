#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void pop_fission_kernel(unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_scanned, spatial_data* fission_particles, float* fission_energy){

	// get tid
	int tid = threadIdx.x+blockIdx.x*blockDim.x;

	// declare shared variables
	__shared__ 	unsigned			n_isotopes;				
	//__shared__ 	unsigned			energy_grid_len;		
	__shared__ 	unsigned			total_reaction_channels;
	__shared__ 	float*				energy_grid;
	__shared__ 	dist_container*		dist_scatter;			
	__shared__ 	dist_container*		dist_energy; 
	__shared__	spatial_data*		space;	
	__shared__	float*				E;		
	__shared__	unsigned*			rn_bank;
	__shared__	unsigned*			yield;	
	__shared__	unsigned*			index;	
	__shared__	unsigned*			isonum;	

	// have thread 0 of block copy all pointers and static info into shared memory
	if (threadIdx.x == 0){
		n_isotopes					= d_xsdata[0].n_isotopes;								
		//energy_grid_len				= d_xsdata[0].energy_grid_len;				
		total_reaction_channels		= d_xsdata[0].total_reaction_channels;
		energy_grid 				= d_xsdata[0].energy_grid;
		dist_scatter 				= d_xsdata[0].dist_scatter;
		dist_energy 				= d_xsdata[0].dist_energy; 
		space						= d_particles[0].space;
		E							= d_particles[0].E;
		rn_bank						= d_particles[0].rn_bank;
		yield						= d_particles[0].yield;
		index						= d_particles[0].index;
		isonum						= d_particles[0].isonum;
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
	unsigned	this_tope		=	isonum[  tid];

	// get array position from prefix scan
	unsigned		position		=	d_scanned[tid];

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
		printf("null pointer, energy array in pop_fission!,tid %u\n",tid);
		return;
	}
	// check S data pointers
	if(dist_scatter == 0x0){
		printf("null pointer, scatter array in pop_fission!,tid %u\n",tid);
		return;
	}

	// check second level pointers
	if(dist_scatter[this_dex].lower == 0x0){printf("pop_fission: null pointer dist_scatter.lower! this_dex %u this_E %6.4E tope %u yield %u\n",this_dex,this_E,this_tope,this_yield); return;}
	if(dist_scatter[this_dex].upper == 0x0){printf("pop_fission: null pointer dist_scatter.upper! this_dex %u this_E %6.4E tope %u yield %u\n",this_dex,this_E,this_tope,this_yield); return;}
	if(dist_energy[ this_dex].upper == 0x0){printf("pop_fission: null pointer dist_energy.upper!  this_dex %u this_E %6.4E tope %u yield %u\n",this_dex,this_E,this_tope,this_yield); return;}
	if(dist_energy[ this_dex].lower == 0x0){printf("pop_fission: null pointer dist_energy.lower!  this_dex %u this_E %6.4E tope %u yield %u\n",this_dex,this_E,this_tope,this_yield); return;}

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
	float 		phi, mu, E0, f, rn1, rn2;
	unsigned	this_law, this_len, this_intt, upper_len, lower_len, pre_index, pre_position;
	float		*this_var, *this_cdf, *this_pdf, *upper_var, *lower_var;

	unsigned 	n_columns	= n_isotopes + total_reaction_channels;
	unsigned	this_col	= this_dex % n_columns;
	unsigned	this_row	= (this_dex-this_col) / n_columns;
	float 		E_of_index0	= energy_grid[this_row];
	float		E_of_index1	= energy_grid[this_row+1];

	if(this_E < E_of_index0 | this_E > E_of_index1){printf("energy outside of distributions in pop!!!! this %6.4E row %6.4E row+1 %6.4E \n",this_E,E_of_index0,E_of_index1);}

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
				interpolate_linear_energy( this_E, e0, e1, nu_t0, nu_t1 ) ;
	if( (this_E > e1 | this_E < e0) & (e0 != e1) ){printf("OUTSIDE bounds in pop_fission!   this_E %6.4E e0 %6.4E e1 %6.4E col %u row %u\n",this_E,e0,e1,this_col,this_row);}

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

		}

		// check for null again
		if( this_var == 0x0){printf("!-!-! null pointer  this_var!\n"); return;}
		if( this_cdf == 0x0){printf("!-!-! null pointer  this_cdf!\n"); return;}
		if( this_pdf == 0x0){printf("!-!-! null pointer  this_pdf!\n"); return;}
		if(upper_var == 0x0){printf("!-!-! null pointer upper_var!\n"); return;}
		if(lower_var == 0x0){printf("!-!-! null pointer lower_var!\n"); return;}
			
		// sample dist, passing the parameters/pointers of the sampled delayed/prompt emission data
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
		else if ( this_law == 7 ){   // maxwell spectrum
	
			// get tabulated temperature
			float t0 = edist_lower.var[0];
			float t1 = edist_upper.var[0];
			float U  = edist_lower.cdf[0];
			float e0 = edist_lower.erg;
			float e1 = edist_upper.erg;
			float  T = 0;
			sampled_E = 99999.0;
	
			// interpolate T
			if (e1==e0 | edist_lower.intt==1){  // in top bin, both values are the same
				T = t0;
			}
			else if (edist_lower.intt==2){// lin-lin interpolation
				T  = (t1 - t0)/(e1 - e0) * (this_E - e0) + t0;
			}
			else{
				printf("dont know what to do!\n");
			}
		
			// restriction
			while (sampled_E > this_E - U){
	
				// rejection sample
				rn1 = get_rand(&rn);
				rn2 = get_rand(&rn);
				while ( rn1*rn1+rn2*rn2 > 1.0 ) {
					rn1 = get_rand(&rn);
					rn2 = get_rand(&rn);
				}
			
				// mcnp5 volIII pg 2-43
				sampled_E = -T * (    rn1*rn1*logf(get_rand(&rn)) / (rn1*rn1+rn2*rn2) + logf(get_rand(&rn))   );
	
			}
	
			// isotropic mu/phi
			mu  = 2.0*get_rand(&rn)-1.0;
			phi = 2.0*pi*get_rand(&rn);
	
		}
		else if ( this_law == 9 ){   //evaopration spectrum
	
			// get tabulated temperature
			float t0 = edist_lower.var[0];
			float t1 = edist_upper.var[0];
			float U  = edist_lower.cdf[0];
			float e0 = edist_lower.erg;
			float e1 = edist_upper.erg;
			float  T = 0.0;
			float  m = 0.0;
		
			// interpolate T
			if (e1==e0 | edist_lower.intt==1){  // in top bin, both values are the same
				T = t0;
			}
			else if (edist_lower.intt==2){// lin-lin interpolation
				T  = (t1 - t0)/(e1 - e0) * (this_E - e0) + t0;
			}
			else{
				printf("dont know what to do!\n");
			}
		
			// rejection sample
			m  = (this_E - U)/T;
			e0 = 1.0-expf(-m);
			float x  = -logf(1.0-e0*get_rand(&rn)) - logf(1.0-e0*get_rand(&rn));
			while (  x>m ) {
				x  = -logf(1.0-e0*get_rand(&rn)) - logf(1.0-e0*get_rand(&rn));
			}
		
			// mcnp5 volIII pg 2-43
			sampled_E = T * x;
		
			// isotropic mu/phi
			mu  = 2.0*get_rand(&rn)-1.0;
			phi = 2.0*pi*get_rand(&rn);
	
		}
		else if( this_law == 11 ){  // energy-dependent maxwellian

			// get tabulated parameters
			float a0 = edist_lower.var[0];
			float a1 = edist_upper.var[0];
			float b0 = edist_lower.cdf[0];
			float b1 = edist_upper.cdf[0];
			float U  = edist_lower.pdf[0];
			float e0 = edist_lower.erg;
			float e1 = edist_upper.erg;
			float  a = 0.0;
			float  b = 0.0;
			float  g = 0.0;
			float  c = 0.0;
			sampled_E = 99999.0;
		
			// interpolate T
			if (e1==e0 | edist_lower.intt==1){  // in top bin, both values are the same
				a = a0;
				b = b0;
			}
			else if (edist_lower.intt==2){// lin-lin interpolation
				a  = (a1 - a0)/(e1 - e0) * (this_E - e0) + a0;
				b  = (b1 - b0)/(e1 - e0) * (this_E - e0) + b0;
				c  = 1.0 + a*b/8.0;
				g  = sqrtf( c*c - 1.0 ) + c;
			}
			else{
				printf("dont know what to do!\n");
			}

			// restriction
			while (sampled_E > this_E - U){
	
				// rejection sample
				rn1 = get_rand(&rn);
				rn2 = get_rand(&rn);
				sampled_E = -a*g*logf(rn1);
				c = (1.0-g)*(1.0-logf(rn1)) - logf(rn2);
				while ( c*c > b*sampled_E ) {
					rn1 = get_rand(&rn);
					rn2 = get_rand(&rn);
					sampled_E = -a*g*logf(rn1);
					c = (1.0-g)*(1.0-logf(rn1)) - logf(rn2);
				}
				
			}

			// isotropic mu/phi
			mu  = 2.0*get_rand(&rn)-1.0;
			phi = 2.0*pi*get_rand(&rn);

		}
		else{
			printf("LAW %u NOT HANDLED IN FISSION POP!\n",this_law);
		}

		// check temp array
		if(fission_energy[   data_dex] 				!= 0.0){printf("NONZERO fission_energy[   data_dex]            = % 6.4E \n",fission_energy[   data_dex] 			);}
		if(fission_particles[data_dex].x			!= 0.0){printf("NONZERO fission_particles[data_dex].x          = % 6.4E \n",fission_particles[data_dex].x			);}
		if(fission_particles[data_dex].y			!= 0.0){printf("NONZERO fission_particles[data_dex].y          = % 6.4E \n",fission_particles[data_dex].y			);}
		if(fission_particles[data_dex].z			!= 0.0){printf("NONZERO fission_particles[data_dex].z          = % 6.4E \n",fission_particles[data_dex].z			);}
		if(fission_particles[data_dex].xhat			!= 0.0){printf("NONZERO fission_particles[data_dex].xhat       = % 6.4E \n",fission_particles[data_dex].xhat		);}
		if(fission_particles[data_dex].yhat			!= 0.0){printf("NONZERO fission_particles[data_dex].yhat       = % 6.4E \n",fission_particles[data_dex].yhat		);}
		if(fission_particles[data_dex].zhat			!= 0.0){printf("NONZERO fission_particles[data_dex].zhat       = % 6.4E \n",fission_particles[data_dex].zhat		);}
		if(fission_particles[data_dex].surf_dist	!= 0.0){printf("NONZERO fission_particles[data_dex].surf_dist  = % 6.4E \n",fission_particles[data_dex].surf_dist	);}
		if(fission_particles[data_dex].enforce_BC	!= 0  ){printf("NONZERO fission_particles[data_dex].enforce_BC = %u     \n",fission_particles[data_dex].enforce_BC	);}
		if(fission_particles[data_dex].norm[0]		!= 0.0){printf("NONZERO fission_particles[data_dex].norm[0]    = % 6.4E \n",fission_particles[data_dex].norm[0]		);}
		if(fission_particles[data_dex].norm[1]		!= 0.0){printf("NONZERO fission_particles[data_dex].norm[1]    = % 6.4E \n",fission_particles[data_dex].norm[1]		);}
		if(fission_particles[data_dex].norm[2]		!= 0.0){printf("NONZERO fission_particles[data_dex].norm[2]    = % 6.4E \n",fission_particles[data_dex].norm[2]		);}

		// set data in temp array since GRID-WISE threadsync cannot be done (easily?)!
		fission_energy[     data_dex ] 					= sampled_E;
		fission_particles[  data_dex ].x				= this_x;
		fission_particles[  data_dex ].y				= this_y;
		fission_particles[  data_dex ].z				= this_z;
		fission_particles[  data_dex ].xhat				= sqrtf(1.0-(mu*mu))*cosf(phi);
		fission_particles[  data_dex ].yhat				= sqrtf(1.0-(mu*mu))*sinf(phi); 
		fission_particles[  data_dex ].zhat				= mu;
		fission_particles[  data_dex ].enforce_BC		= 0;
		fission_particles[  data_dex ].surf_dist		= 999999.0;
		
		//if(data_dex<=9){printf("array index %u, E = % 6.4E d_fissile_energy[ data_dex ] = % 6.4E\n",data_dex,sampled_E,E[ data_dex ]);}

	}

	// write current seed out
	rn_bank[tid] = rn;

}

/**
 * \brief a
 * \details b
 *
 * @param[in]    NUM_THREADS        - the number of threads to run per thread block
 * @param[in]    N                  - the total number of threads to launch on the grid
 * @param[in]    d_xsdata           - device pointer to cross section data pointer array
 * @param[in]    d_particles        - device pointer to particle data pointer array
 * @param[in]    d_scanned          - device pointer to array of the cumulative sum (scan) of the yield array, used to find final index where new particles will be written
 * @param[in]    fission_particles  - device pointer to intermadiate spatial data array where popped values will be written
 * @param[in]    fission_energy     - device pointer to intermadiate energy data array where popped values will be written
 */ 
void pop_fission( unsigned NUM_THREADS, unsigned N, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_scanned, spatial_data* fission_particles, float* fission_energy ){

	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	pop_fission_kernel <<< blks, NUM_THREADS >>> ( N, d_xsdata, d_particles, d_scanned, fission_particles, fission_energy);
	check_cuda(cudaThreadSynchronize());

}

