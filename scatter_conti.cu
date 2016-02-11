#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "warp_device.cuh"
#include "check_cuda.h"


__global__ void scatter_conti_kernel(unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

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

	// return immediately if out of bounds
	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	//remap to active
	int tid				=	d_remap[starting_index + tid_in];
	unsigned this_rxn 	=	rxn[    starting_index + tid_in];

	// print and return if wrong
	if ( this_rxn!=91 ){printf("level scattering kernel accessing wrong reaction @ dex %u rxn %u\n",tid, this_rxn);return;} 

	// check E data pointers
	if(dist_energy == 0x0){
		printf("null pointer, energy array in continuum scatter!,tid %u rxn %u\n",tid,this_rxn);
		return;
	}

	//constants
	//const float  	pi			=   3.14159265359;
	const float  	m_n			=   1.00866491600;		// u
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
	float 		mu, E0, A, R;
	unsigned 	dist_index[1];    // must be declared this way in order to write to passed pointer, why??

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
		sampled_E = scale_to_bins(	f, E0, 
									 this_edist.var[0],  this_edist.var[ this_edist.len-1], 
									edist_lower.var[0], edist_lower.var[edist_lower.len-1], 
									edist_upper.var[0], edist_upper.var[edist_upper.len-1] );

		// sample mu isotropically
		mu  = 2.0*get_rand(&rn)-1.0;

	}
//	else if (law==9){   //evaopration spectrum
//
//		// get tabulated temperature
//		float t0 = this_Earray[ offset              ];
//		float t1 = this_Earray[ offset + 1          ];
//		float U  = this_Earray[ offset + vlen       ];
//			  e0 = this_Earray[ offset + vlen*2     ];
//			  e1 = this_Earray[ offset + vlen*2 + 1 ];
//		float  T = 0.0;
//		float  m = 0.0;
//
//		// interpolate T
//		if (e1==e0){  // in top bin, both values are the same
//			T = t0;
//		}
//		else if (intt==2){// lin-lin interpolation
//			m = (this_E - e0)/(e1 - e0);
//			T = (1.0 - m)*t0 + m*t1;
//		}
//		else if(intt==1){// histogram interpolation
//			T  = (t1 - t0)/(e1 - e0) * this_E + t0;
//		}
//
//		// rejection sample
//		m  = (this_E - U)/T;
//		e0 = 1.0-expf(-m);
//		float x  = -logf(1.0-e0*get_rand(&rn)) - logf(1.0-e0*get_rand(&rn));
//		while (  x>m ) {
//			x  = -logf(1.0-e0*get_rand(&rn)) - logf(1.0-e0*get_rand(&rn));
//		}
//
//		// mcnp5 volIII pg 2-43
//		sampled_E = T * x;
//
//		//isotropic mu
//		if (this_Sarray==0x0){
//			mu  = 2.0*get_rand(&rn)-1.0;
//		}
//		else{
//			printf("law 9 in cscatter has angular tables\n");
//		}
//
//	}
	else if (this_law==44){

		// make sure scatter array is present
		if(dist_scatter == 0x0){
			printf("null pointer, scatter array in continuum !,dex %u rxn %u tope %u E %6.4E \n",this_dex,this_rxn,this_tope,this_E);
			return;
		}

		// correct if below lower energy?
		//if(this_E<this_edist.var[this_edist.len-1]){
		//	printf("above last e! \n");
		//	//this_E = last_E;
		//}

		// compute interpolation factor
		if(f<0){
			printf("DATA NOT WITHIN ENERGY INTERVAL tid %u rxn %u\n",tid,this_rxn);
		}

		// sample tabular on energy, but get index as well as value
		E0 = sample_continuous_tablular( 	dist_index ,
											this_edist.len , 
											2 , 
											get_rand(&rn) , 
											this_edist.var , 
											this_edist.pdf, 
											this_edist.cdf );
		//scale it to bins 
		sampled_E = scale_to_bins(	f, E0, 
									 this_edist.var[0],  this_edist.var[ this_edist.len-1], 
									edist_lower.var[0], edist_lower.var[edist_lower.len-1], 
									edist_upper.var[0], edist_upper.var[edist_upper.len-1] );

		// find correlated mu
		if (this_sdist.intt==1){
			A	=	this_sdist.var[dist_index[0]];
			R	=	this_sdist.cdf[dist_index[0]];
		}
		else if (this_sdist.intt==2){
			A	=	interpolate_linear_energy(	E0,
												this_edist.var[dist_index[0]],
												this_edist.var[dist_index[0]+1],
												this_sdist.var[dist_index[0]],
												this_sdist.var[dist_index[0]+1]);
			R	=	interpolate_linear_energy(	E0,
												this_edist.var[dist_index[0]],
												this_edist.var[dist_index[0]+1],
												this_sdist.cdf[dist_index[0]],
												this_sdist.cdf[dist_index[0]+1]);
		}
		else{
			printf("INTT=%u NOT HANDLED in law %u of continuum scatter!",this_sdist.law,this_sdist.intt);
		}
		float rn1 	= get_rand(&rn);
		if( get_rand(&rn)>R ){
			float T = (2.0*rn1-1.0)*sinhf(A);
			mu		= logf(T+sqrtf(T*T+1.0))/A;
		}
		else{
			mu		= logf(rn1*expf(A)+(1.0-rn1)*expf(-A))/A;
		}

	}
//	else if (law==61){
//
//		unsigned distloc, vloc;
//		float r = (this_E-last_E)/(next_E-last_E);
//		last_e_start = this_Earray[ offset ];
//		last_e_end   = this_Earray[ offset + vlen - 1 ];
//		next_e_start = this_Earray[ offset + 3*vlen ];
//		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
//	
//		rn1 = get_rand(&rn);
//		rn2 = get_rand(&rn);
//	
//		//sample energy dist
//		sampled_E = 0.0;
//		if(  rn2 >= r ){   //sample last E
//			distloc = 1;   // use the first flattened array
//			diff = next_e_end - next_e_start;
//			e_start = next_e_start;
//			for ( n=0 ; n<vlen-1 ; n++ ){
//				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
//				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
//				pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
//				pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
//				e0  		= this_Earray[ (offset          ) + n+0];
//				e1  		= this_Earray[ (offset          ) + n+1]; 
//				if( rn1 >= cdf0 & rn1 < cdf1 ){
//					break;
//				}
//			}
//		}
//		else{
//			distloc = this_Sarray[0];   // get location of the next flattened array
//			diff = next_e_end - next_e_start;
//			e_start = next_e_start;
//			for ( n=0 ; n<next_vlen-1 ; n++ ){
//				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
//				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
//				pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
//				pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
//				e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
//				e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
//				if( rn1 >= cdf0 & rn1 < cdf1 ){
//					break;
//				}
//			}
//		}
//	
//		if (intt==2){// lin-lin interpolation
//			float m 	= (pdf1 - pdf0)/(e1-e0);
//			float arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
//			if(arg<0){
//				E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
//			}
//			else{
//				E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
//			}
//		}
//		else if(intt==1){// histogram interpolation
//			E0 = e0 + (rn1-cdf0)/pdf0;
//		}
//		
//		//scale it
//		E1 = last_e_start + r*( next_e_start - last_e_start );
//		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
//		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
//
//		//
//		// sample mu from tabular distributions
//		//
//
//		// get parameters
//		unsigned vlen_S ;
//		if(distloc){
//			unsigned l = this_Sarray[0];
//			vloc   = this_Sarray[l + n] + (l + next_vlen) ; // get appropriate vector location for this E_out
//			}                
//		else{   
//			vloc   = this_Sarray[1 + n] + (1 + vlen) ;     
//		}
//		vlen_S = this_Sarray[vloc + 0];        // vector length
//		intt   = this_Sarray[vloc + 1];        // interpolation type
//		//printf("distloc %u vloc %u vlen_S %u intt %u \n",distloc,vloc,vlen_S,intt);
//
//		// sample the dist
//		rn1 = get_rand(&rn);
//		for ( n=0 ; n<vlen-1 ; n++ ){
//			cdf0 		= this_Sarray[ (vloc + 2 +   vlen_S ) + n+0];
//			cdf1  		= this_Sarray[ (vloc + 2 +   vlen_S ) + n+1];
//			pdf0		= this_Sarray[ (vloc + 2 + 2*vlen_S ) + n+0];
//			pdf1		= this_Sarray[ (vloc + 2 + 2*vlen_S ) + n+1];
//			e0   		= this_Sarray[ (vloc + 2            ) + n+0];
//			e1   		= this_Sarray[ (vloc + 2            ) + n+1];
//			if( rn1 >= cdf0 & rn1 < cdf1 ){
//				break;
//			}
//		}
//
//		// interpolate
//		if (e1==e0){  // in top bin, both values are the same
//				mu = e0;
//			}
//		else if (intt==2){// lin-lin interpolation
//			r = (rn1 - cdf0)/(cdf1 - cdf0);
//           mu = (1.0 - r)*e0 + r*e1;
//		}
//		else if(intt==1){// histogram interpolation
//			mu  = (e1 - e0)/(cdf1 - cdf0) * rn1 + e0;
//		}
//		else{
//			printf("intt in law 61 in cscatter is invlaid (%u)!\n",intt);
//		}
//		
//
//	}
	else{

		printf("LAW %u NOT HANDLED IN CONTINUUM SCATTER!  rxn %u\n",this_law,this_rxn);

	}

	// check errors
	if (!isfinite(sampled_E) | sampled_E < 0.0){
		printf("continuum scatter mis-sampled tid %i data_dex %u E %6.4E dist len %u dist_index %u... \n",tid_in,tid,sampled_E,this_sdist.len,dist_index[0]);
	}
	if (!isfinite(mu) | mu < -1.0 | mu > 1.0){
		printf("continuum scatter mis-sampled tid %i data_dex %u mu %6.4E dist len %u dist_index %u... \n",tid_in,tid,mu,this_sdist.len,dist_index[0]);
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

	//printf("tid %d law %u sampled_E %6.4E mu %6.4E\n",tid,this_law,sampled_E,mu);
	
	// write universal results
	E[tid]			=	E_new;
	space[tid].xhat	=	hats_new.x;
	space[tid].yhat	=	hats_new.y;
	space[tid].zhat	=	hats_new.z;
	rn_bank[tid]	=	rn;

}

void scatter_conti( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index, cross_section_data* d_xsdata, particle_data* d_particles, unsigned* d_remap){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;
	
	scatter_conti_kernel <<< blks, NUM_THREADS , 0 , stream >>> ( N, starting_index, d_xsdata, d_particles, d_remap );
	check_cuda(cudaThreadSynchronize());

}

