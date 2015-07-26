#include <cuda.h>
#include <stdio.h>
#include "datadef.h"
#include "wfloat3.h"
#include "LCRNG.cuh"
#include "binary_search.h"

__device__ void process_secondaries(source_point* space_out, float* E_out, unsigned* rn, unsigned this_tope, float this_awr, float* this_Earray, float* this_Sarray){

	// internal data
	unsigned 	n, offset, vlen, next_vlen, law, intt;
	float 		sampled_E, phi, mu, rn1, rn2, last_E, next_E, e_start, E0, E1, Ek, next_e_start, next_e_end, last_e_start, last_e_end, diff, r;
	float 		cdf0, cdf1, e0, e1, m, pdf0, pdf1, arg, x,y,z, A, R;
	const float pi 		= 3.14159265359;
	const float Emin 	= 1e-11;
	const float Emax 	= 200.0;
	float this_E = E_out[0];

	//read in values
	offset = 6;
	memcpy(&last_E,   	&this_Earray[0], sizeof(float));
	memcpy(&next_E,   	&this_Earray[1], sizeof(float));
	memcpy(&vlen,   	&this_Earray[2], sizeof(float));
	memcpy(&next_vlen,	&this_Earray[3], sizeof(float));
	memcpy(&law, 		&this_Earray[4], sizeof(float)); 
	memcpy(&intt, 		&this_Earray[5], sizeof(float)); 


	if (law == 4 | law == 44 | law == 61) {
		r = (this_E-last_E)/(next_E-last_E);
		last_e_start = this_Earray[ offset ];
		last_e_end   = this_Earray[ offset + vlen - 1 ];
		next_e_start = this_Earray[ offset + 3*vlen ];
		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
	}
	
	if (law==4){ // tabular
		rn1 = get_rand(rn);
		rn2 = get_rand(rn);
	
		//sample energy dist
		sampled_E = 0.0;
		if(  rn2 >= r ){   //sample last E
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			for ( n=0 ; n<vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
				e0  		= this_Earray[ (offset          ) + n+0];
				e1  		= this_Earray[ (offset          ) + n+1]; 
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
		}
		else{
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			for ( n=0 ; n<next_vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
				e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
				e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
		}
	
		if (intt==2){// lin-lin interpolation
			m 	= (pdf1 - pdf0)/(e1-e0);
			arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
			if(arg<0){
				E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
			}
			else{
				E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
			}
		}
		else if(intt==1){// histogram interpolation
			E0 = e0 + (rn1-cdf0)/pdf0;
		}
		
		//scale it
		E1 = last_e_start + r*( next_e_start - last_e_start );
		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;

		//isotropic mu
		mu  = 2.0*get_rand(rn)-1.0;

	}
	else if (law==44){
		if(this_Sarray == 0x0){  //Sarray value is nu for 918, might have to put check in here
			printf("null pointer in pop multiplicity Sarray!, tope %u E %6.4E\n",this_tope,this_E);
			return;
		}
		//printf("%6.4E %6.4E %u %u %u ... %6.4E %6.4E ... %6.4E %6.4E %6.4E\n",this_Earray[0],this_Earray[1],vlen,next_vlen,law,this_Earray[5],this_Earray[6],this_Earray[offset],this_Earray[offset+vlen],this_Earray[offset+vlen+1]);
		if(this_E<last_E){this_E=last_E;}
		float r = (this_E-last_E)/(next_E-last_E);
		if(r<0){
			printf("r less than zero in source pop for multiplicity, r % 10.8E isotope %u this_E % 10.8E last_E % 10.8E next_E % 10.8E\n",r,this_tope,this_E,last_E,next_E);
		}
		last_e_start = this_Earray[ offset ];
		last_e_end   = this_Earray[ offset + vlen - 1 ];
		next_e_start = this_Earray[ offset + 3*vlen ];
		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
	
		//sample energy dist
		sampled_E = 0.0;
		rn1 = get_rand(rn);
		if(  get_rand(rn) >= r ){   //sample last E
			diff = last_e_end - last_e_start;
			e_start = last_e_start;
			//n = binary_search( &this_Earray[ offset + vlen ] , rn1, vlen);
			for ( n=0 ; n<vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
				e0  		= this_Earray[ (offset          ) + n+0];
				e1  		= this_Earray[ (offset          ) + n+1]; 
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
			offset = 6;
			A = this_Sarray[ (offset)      + n ];
			R = this_Sarray[ (offset+vlen) + n ];
		}
		else{
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			//n = binary_search( &this_Earray[ offset + 3*vlen + next_vlen] , rn1, next_vlen);
			for ( n=0 ; n<next_vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
				e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
				e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
			offset = 6;
			A = this_Sarray[ (offset+3*vlen)           +n  ] ;
			R = this_Sarray[ (offset+3*vlen+next_vlen) +n  ];
		}
		
		// histogram interpolation, intt=1
		E0 = e0 + (rn1-cdf0)/pdf0;
		
		//scale it
		E1 = last_e_start + r*( next_e_start - last_e_start );
		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
		
		// find mu
		rn1 = get_rand(rn);
		if(get_rand(rn)>R){
			float T = (2.0*rn1-1.0)*sinhf(A);
			mu = logf(T+sqrtf(T*T+1.0))/A;
		}
		else{
				mu = logf(rn1*expf(A)+(1.0-rn1)*expf(-A))/A;
			}
	}
	else if (law==61){

		unsigned distloc, vloc;
		float r = (this_E-last_E)/(next_E-last_E);
		last_e_start = this_Earray[ offset ];
		last_e_end   = this_Earray[ offset + vlen - 1 ];
		next_e_start = this_Earray[ offset + 3*vlen ];
		next_e_end   = this_Earray[ offset + 3*vlen + next_vlen - 1];
	
		rn1 = get_rand(rn);
		rn2 = get_rand(rn);
	
		//sample energy dist
		sampled_E = 0.0;
		if(  rn2 >= r ){   //sample last E
			distloc = 1;   // use the first flattened array
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			for ( n=0 ; n<vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset +   vlen ) + n+0];
				cdf1 		= this_Earray[ (offset +   vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 2*vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 2*vlen ) + n+1];
				e0  		= this_Earray[ (offset          ) + n+0];
				e1  		= this_Earray[ (offset          ) + n+1]; 
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
		}
		else{
			distloc = this_Sarray[0];   // get location of the next flattened array
			diff = next_e_end - next_e_start;
			e_start = next_e_start;
			for ( n=0 ; n<next_vlen-1 ; n++ ){
				cdf0 		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+0];
				cdf1  		= this_Earray[ (offset + 3*vlen +   next_vlen ) + n+1];
				pdf0		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+0];
				pdf1		= this_Earray[ (offset + 3*vlen + 2*next_vlen ) + n+1];
				e0   		= this_Earray[ (offset + 3*vlen               ) + n+0];
				e1   		= this_Earray[ (offset + 3*vlen               ) + n+1];
				if( rn1 >= cdf0 & rn1 < cdf1 ){
					break;
				}
			}
		}
	
		if (intt==2){// lin-lin interpolation
			float m 	= (pdf1 - pdf0)/(e1-e0);
			float arg = pdf0*pdf0 + 2.0 * m * (rn1-cdf0);
			if(arg<0){
				E0 = e0 + (e1-e0)/(cdf1-cdf0)*(rn1-cdf0);
			}
			else{
				E0 	= e0 + (  sqrtf( arg ) - pdf0) / m ;
			}
		}
		else if(intt==1){// histogram interpolation
			E0 = e0 + (rn1-cdf0)/pdf0;
		}
		
		//scale it
		E1 = last_e_start + r*( next_e_start - last_e_start );
		Ek = last_e_end   + r*( next_e_end   - last_e_end   );
		sampled_E = E1 +(E0-e_start)*(Ek-E1)/diff;
	
		//
		// sample mu from tabular distributions
		//
	
		// get parameters
		unsigned vlen_S ;
		if(distloc){
			unsigned l = this_Sarray[0];
			vloc   = this_Sarray[l + n] + (l + next_vlen) ; // get appropriate vector location for this E_out
			}                
		else{   
			vloc   = this_Sarray[1 + n] + (1 + vlen) ;     
		}
		vlen_S = this_Sarray[vloc + 0];        // vector length
		intt   = this_Sarray[vloc + 1];        // interpolation type
		//printf("distloc %u vloc %u vlen_S %u intt %u \n",distloc,vloc,vlen_S,intt);
	
		// sample the dist
		rn1 = get_rand(rn);
		for ( n=0 ; n<vlen-1 ; n++ ){
			cdf0 		= this_Sarray[ (vloc + 2 +   vlen_S ) + n+0];
			cdf1  		= this_Sarray[ (vloc + 2 +   vlen_S ) + n+1];
			pdf0		= this_Sarray[ (vloc + 2 + 2*vlen_S ) + n+0];
			pdf1		= this_Sarray[ (vloc + 2 + 2*vlen_S ) + n+1];
			e0   		= this_Sarray[ (vloc + 2            ) + n+0];
			e1   		= this_Sarray[ (vloc + 2            ) + n+1];
			if( rn1 >= cdf0 & rn1 < cdf1 ){
				break;
			}
		}
	
		// interpolate
		if (e1==e0){  // in top bin, both values are the same
				mu = e0;
			}
		else if (intt==2){// lin-lin interpolation
			r = (rn1 - cdf0)/(cdf1 - cdf0);
 	       mu = (1.0 - r)*e0 + r*e1;
		}
		else if(intt==1){// histogram interpolation
			mu  = (e1 - e0)/(cdf1 - cdf0) * rn1 + e0;
		}
		else{
			printf("intt in law 61 in cscatter is invlaid (%u)!\n",intt);
		}
	
	}
	else if(law==7){   // maxwellian fission

		// get tabulated temperature
		float t0 = this_Earray[ offset     ];
		float t1 = this_Earray[ offset + 1 ];
		float U  = this_Earray[ offset + vlen       ];
		      e0 = this_Earray[ offset + vlen*2     ];
		      e1 = this_Earray[ offset + vlen*2 + 1 ];
		float  T = 0.0;

		// interpolate T
		if (e1==e0){  // in top bin, both values are the same
			T = t0;
		}
		else if (intt==2){// lin-lin interpolation
			m = (this_E - e0)/(e1 - e0);
            T = (1.0 - m)*t0 + m*t1;
		}
		else if(intt==1){// histogram interpolation
			T  = (t1 - t0)/(e1 - e0) * this_E + t0;
		}
		else{
			printf("something is wrong in pop for law 7");
		}

		// rejection sample
		sampled_E = 9999999.9;
		while ( sampled_E > (this_E - U)){
			m = cosf(pi*get_rand(rn)/2.0);
			sampled_E = -T * ( m*m*logf(get_rand(rn))  +   logf(get_rand(rn)) );
		}

		//isotropic mu
		mu  = 2.0*get_rand(rn)-1.0;

	}
	else if (law==9){   //evaporation spectrum

		// get tabulated temperature  
		float t0 = this_Earray[ offset              ];
		float t1 = this_Earray[ offset + 1          ];
		float U  = this_Earray[ offset + vlen       ];
		      e0 = this_Earray[ offset + vlen*2     ];
		      e1 = this_Earray[ offset + vlen*2 + 1 ];
		float  T = 0.0;

		// interpolate T
		if (e1==e0){  // in top bin, both values are the same
			T = t0;
		}
		else if (intt==2){// lin-lin interpolation
			m = (this_E - e0)/(e1 - e0);
            T = (1.0 - m)*t0 + m*t1;
		}
		else if(intt==1){// histogram interpolation
			T  = (t1 - t0)/(e1 - e0) * this_E + t0;
		}

		// rejection sample
		m  = (this_E - U)/T;
		e0 = 1.0-expf(-m);
		x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
		while (  x>m ) {
			x  = -logf(1.0-e0*get_rand(rn)) - logf(1.0-e0*get_rand(rn));
		}

		// mcnp5 volIII pg 2-43
		sampled_E = T * x;

		//isotropic mu
		mu  = 2.0*get_rand(rn)-1.0;

	}
	else if(law==66){   // N-body phase space

		// get tabulated temperature
		float Q 		= this_Earray[ offset     ];
		float nbodies 	= this_Earray[ offset + 1 ];
		if (nbodies>3) {printf("nbodies in law 66  is greater than 3!!!\n");}
		float A 		= this_Earray[ offset + 2 ];
		unsigned this_yield1 = 2; //ONLY FOR 1002 n,2n!
		float Emax 		= (this_yield1-1.0)/this_yield1 * (A/(A+1)*this_E + Q); 
		float rn1       = get_rand(rn);
		float rn2       = get_rand(rn);
		float rn3       = get_rand(rn);
		float rn4       = get_rand(rn);

		//// rejection sample
		while ( (rn1*rn1 + rn2*rn2) > 1.0 ){
			rn1 = get_rand(rn);
			rn2 = get_rand(rn);
		}
		while ( (rn3*rn3 + rn4*rn4) > 1.0 ){
			rn3 = get_rand(rn);
			rn4 = get_rand(rn);
		}

		float x = -rn1 * logf(rn1*rn1 + rn2*rn2)/ (rn1*rn1 + rn2*rn2) - logf(get_rand(rn)) ;
		float y = -rn3 * logf(rn3*rn3 + rn4*rn4)/ (rn3*rn3 + rn4*rn4) - logf(get_rand(rn)) ;

		sampled_E = Emax * x /( x + y );

		//isotropic mu
		mu  = 2.0*get_rand(rn)-1.0;

	}
	else{
		printf("LAW %u NOT HANDLED IN FISSION!\n",law);
	}

	//sample isotropic directions
	rn1 = get_rand(rn);
	rn2 = get_rand(rn);
	mu  = 2.0*rn1-1.0; 
	phi = 2.0*pi*rn2;
	x = sqrtf(1.0-(mu*mu))*cosf(phi);
	y = sqrtf(1.0-(mu*mu))*sinf(phi);
	z = mu;
	
	//check limits
	if (sampled_E >= Emax){sampled_E = Emax * 0.99;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}
	if (sampled_E <= Emin){sampled_E = Emin * 1.01;}//printf("enforcing limits in pop data_dex=%u, sampled_E = %6.4E\n",data_dex,sampled_E);}

	// set data
	space_out[0].xhat 		= x;
	space_out[0].yhat 		= y;
	space_out[0].zhat 		= z;
	space_out[0].enforce_BC = 0;
	space_out[0].surf_dist 	= 99999.0;
	space_out[0].macro_t 	= 8.675309;
	E_out 	 [0] 			= sampled_E;

}

__global__ void fission_kernel(unsigned N, unsigned starting_index , unsigned* remap, unsigned* isonum, unsigned* index, unsigned* rn_bank, float* E, source_point* space, unsigned* rxn, float* awr_list, unsigned* yield, float* weight, float** scatterdata, float** energydata){

	int tid_in = threadIdx.x+blockIdx.x*blockDim.x;
	if (tid_in >= N){return;}

	//remap to active
	int      tid      = remap[starting_index + tid_in];
	unsigned this_rxn = rxn  [starting_index + tid_in];

	// check
	unsigned 		this_yield 	= yield[tid];
	if (this_yield!=0){
		printf("tid %u REACTION %u IS IN FISSION WITH NONZERO YIELD %u!\n",tid,this_rxn,this_yield);
		return;
	}

	// load in more external data
	unsigned 		this_tope 	= isonum [tid];
	unsigned 		dex  		= index  [tid];
	unsigned 		rn 			= rn_bank[tid];
	float 			this_weight = weight [tid]; 
	float*	 		this_Sarray = scatterdata[dex];
	float*			this_Earray = energydata [dex];
	unsigned 		inu 		= 0;
	float 			nu 			= 0.0;

	// print and return if wrong
	if (this_rxn < 811 | this_rxn > 845){
		printf("fission kernel accessing wrong reaction @ remapped dex %u sorted dex %u rxn %u\n",tid,starting_index + tid_in, this_rxn);
		return;
	} 

	// determine integer yields for fission
	if (this_rxn == 818 | this_rxn == 819 | this_rxn == 820 | this_rxn == 821){
	
		//load nu value, since e search has alrady been done, nu should be where the scatter array is (fission is always isotropic)
		memcpy(&nu, &scatterdata[dex], sizeof(float));
		if (nu==0.0){
			nu=2.8;
			printf("something is wrong with fission yields, nu = %6.4E, guessing %4.2f\n",0.0,nu); 
		}

		//  multiply nu by weight
		nu = this_weight * nu;

		// get integer part
		inu = (unsigned) nu;
		
		// sample floor or ceil based on fractional part
		if((float)inu+get_rand(&rn) <= nu){
			this_yield = inu+1;
		}
		else{
			this_yield = inu;
		}

		// put in 900 to terminate
 		this_rxn += 100;

	}
	else{  // needs to be resampled as scatter, adjust weight to yield

		// check data array pointers
		if(this_Earray == 0x0){
			printf("null pointer in multiplicity Earray,tid %u dex %u rxn %u tope %u E %6.4E\n",tid,dex,this_rxn,this_tope,E[tid]);
			return;
		}

		// fission yield is zero
		this_yield = 0;

		// set weights as multiplicity
		if(this_rxn == 817 | this_rxn == 825 | this_rxn == 842){
			this_weight = 3.0*this_weight;  
		}
		else if(this_rxn == 816 | this_rxn==824 | this_rxn == 811 | this_rxn == 824 | this_rxn == 829 | this_rxn == 830 | this_rxn == 841){
			this_weight = 2.0*this_weight;  
		}
		else{
			this_weight = 1.0*this_weight;
		}

		// determine new E/dir, values written internally
		process_secondaries(&space[tid], &E[tid], &rn, this_tope, awr_list[this_tope], this_Earray, this_Sarray);

	}

	if(this_rxn==818){printf("rxn still 818 somehow!\n");}

	// write 
	weight[tid] = this_weight;
	yield[tid]  = this_yield;
	rn_bank[tid] = rn;  
	rxn[starting_index + tid_in] = this_rxn;

}


void fission( cudaStream_t stream, unsigned NUM_THREADS, unsigned N, unsigned starting_index , unsigned* d_remap, unsigned* d_isonum, unsigned* d_index, unsigned* d_rn_bank, float* d_E, source_point* d_space, unsigned* d_rxn, float* d_awr_list, unsigned* d_yield, float* d_weight, float** scatterdata, float** energydata){

	if(N<1){return;}
	unsigned blks = ( N + NUM_THREADS - 1 ) / NUM_THREADS;

	fission_kernel <<< blks, NUM_THREADS , 0 , stream >>> ( N, starting_index , d_remap, d_isonum, d_index, d_rn_bank, d_E, d_space, d_rxn, d_awr_list, d_yield, d_weight, scatterdata, energydata);
	cudaThreadSynchronize();

}

